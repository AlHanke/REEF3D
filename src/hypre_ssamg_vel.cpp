/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "hypre_ssamg.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "fieldint4.h"

// ====================================================================================
// Implicit momentum-diffusion (Helmholtz) solves on the staggered velocity systems.
//
// var 1/2/3 -> u/v/w. The velocity fields are stored as CELL-CENTRED MultiFabs (the
// staggering is implicit, see field_amrex), so they share pressure's cell-centred 7-point
// grid -- make_grid_7p is reused, no separate grid is built. The linear system A x = b comes
// straight from a->M / a->rhsvec as written by idiff2_FS* in a single ULOOP/VLOOP/WLOOP pass;
// number_velocity rebuilds that exact enumeration so row n maps back to its cell.
//
// This is a single-level path: the Helmholtz operator is non-singular (CPOR/(alpha*dt) on the
// diagonal makes it diagonally dominant), so unlike pressure there is NO coarse-fine coupling,
// NO nullspace projection and NO reference pin. Multi-level diffusion (a C-F Helmholtz stencil)
// is a separate, future task; guarded here by the nlevs loop doing nothing extra.
// ====================================================================================

void hypre_ssamg::start_solver123(lexer *p, fdm *a, ghostcell *pgc, field &f, int var)
{
    p->solveriter = 0;

    #if USE_AMREX
    if(p->changed)
    {
        destroy_grid();
        make_grid_7p(p, a, pgc);
    }
    #endif

    #if USE_AMREX
    if(!solver_created || created_nlevs != p->nlevs)
    #else
    if(!solver_created)
    #endif
    {
        delete_solver();
        create_solver(p, pgc);
    }

    fill_matrix_vel(p, a, pgc, f, var);

    solve(p);

    fillbackvec_vel(p, f, var);

    std::fill(a->rhsvec.V.begin(), a->rhsvec.V.end(), 0.0);
}

void hypre_ssamg::number_velocity(lexer *p, fieldint4 &cval, int var)
{
    // -1 = un-numbered (identity row); otherwise the sequential index that idiff2_FS* used
    // when it wrote a->M[n]/a->rhsvec[n]. Must use the SAME loop macro (and last-face trim)
    // as the diffusion fill, or the matrix rows map to the wrong cells.
    cval.setVal(-1, true);

    count = 0;
    if(var==1)
    {
        ULOOP { cval(i,j,k) = count; ++count; }
    }
    else if(var==2)
    {
        VLOOP { cval(i,j,k) = count; ++count; }
    }
    else
    {
        WLOOP { cval(i,j,k) = count; ++count; }
    }
}

void hypre_ssamg::fill_matrix_vel(lexer *p, fdm *a, ghostcell *pgc, field &f, int var)
{
    fieldint4 cval(p);
    number_velocity(p, cval, var);

    nentries = stencil_size;
    for (int e = 0; e < nentries; ++e)
        stencil_indices[e] = e;

    // The 2D solvers (idiff2_FSu_2D / idiff2_FSw_2D) are a 5-point stencil and never write
    // the y-couplings M.e/M.w, so those slots hold stale values from a previous a->M use.
    // Force them to zero unless this is a 3D (y-active) run, where idiff2_FS sets them.
    const bool ydir = (p->j_dir != 0);

#if USE_AMREX
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto&       cval_mf  = cval.GetMultiFab(lev);
        const auto& field_mf = f.GetMultiFab(lev);

        for (amrex::MFIter mfi(cval_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto cval_arr  = cval_mf.const_array(mfi);
            const auto field_arr = field_mf.const_array(mfi);

            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            // Matrix A: 7 stencil entries per cell (order p, s(x-), n(x+), e(y-), w(y+), b(z-), t(z+)).
            values.resize(ncells * 7);
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                const int n = cval_arr(ii, jj, kk);
                if (n >= 0)
                {
                    values[cnt++] = a->M.p[n];
                    values[cnt++] = a->M.s[n];
                    values[cnt++] = a->M.n[n];
                    values[cnt++] = ydir ? a->M.e[n] : 0.0;
                    values[cnt++] = ydir ? a->M.w[n] : 0.0;
                    values[cnt++] = a->M.b[n];
                    values[cnt++] = a->M.t[n];
                }
                else
                {
                    values[cnt++] = 1.0;
                    values[cnt++] = 0.0; values[cnt++] = 0.0;
                    values[cnt++] = 0.0; values[cnt++] = 0.0;
                    values[cnt++] = 0.0; values[cnt++] = 0.0;
                }
            }
            HYPRE_SStructMatrixSetBoxValues(A, lev, lo, hi, variable, nentries, stencil_indices, values.data());

            // Initial guess x: the incoming field (idiff seeds it with u_in via CopyFrom).
            values.resize(ncells);
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
                values[cnt++] = (cval_arr(ii, jj, kk) >= 0) ? field_arr(ii, jj, kk) : 0.0;
            HYPRE_SStructVectorSetBoxValues(x, lev, lo, hi, variable, values.data());

            // RHS b
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                const int n = cval_arr(ii, jj, kk);
                values[cnt++] = (n >= 0) ? a->rhsvec.V[n] : 0.0;
            }
            HYPRE_SStructVectorSetBoxValues(b, lev, lo, hi, variable, values.data());
        }
    }
#else
    values.resize(p->knox * p->knoy * p->knoz * 7);

    count = 0;
    KJILOOP
    {
        const int n = cval(i, j, k);
        if (n >= 0)
        {
            values[count++] = a->M.p[n];
            values[count++] = a->M.s[n];
            values[count++] = a->M.n[n];
            values[count++] = ydir ? a->M.e[n] : 0.0;
            values[count++] = ydir ? a->M.w[n] : 0.0;
            values[count++] = a->M.b[n];
            values[count++] = a->M.t[n];
        }
        else
        {
            values[count++] = 1.0;
            values[count++] = 0.0; values[count++] = 0.0;
            values[count++] = 0.0; values[count++] = 0.0;
            values[count++] = 0.0; values[count++] = 0.0;
        }
    }
    HYPRE_SStructMatrixSetBoxValues(A, 0, ilower, iupper, variable, nentries, stencil_indices, values.data());

    count = 0;
    KJILOOP
    {
        const int n = cval(i, j, k);
        values[count++] = (n >= 0) ? f(i, j, k) : 0.0;
    }
    HYPRE_SStructVectorSetBoxValues(x, 0, ilower, iupper, variable, values.data());

    count = 0;
    KJILOOP
    {
        const int n = cval(i, j, k);
        values[count++] = (n >= 0) ? a->rhsvec.V[n] : 0.0;
    }
    HYPRE_SStructVectorSetBoxValues(b, 0, ilower, iupper, variable, values.data());
#endif

    HYPRE_SStructMatrixAssemble(A);
    HYPRE_SStructVectorAssemble(x);
    HYPRE_SStructVectorAssemble(b);

#if USE_AMREX
    // Single-level velocity solves keep object_type==HYPRE_SSTRUCT, so this is a no-op here;
    // kept for symmetry with fill_matrix4 should multi-level diffusion ever extract ParCSR.
    if (p->nlevs > 1)
    {
        HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
        HYPRE_SStructVectorGetObject(b, (void**) &par_b);
        HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    }
#endif
}

void hypre_ssamg::fillbackvec_vel(lexer *p, field &f, int var)
{
    fieldint4 cval(p);
    number_velocity(p, cval, var);

#if USE_AMREX
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto&       field_mf = f.GetMultiFab(lev);
        const auto& cval_mf  = cval.GetMultiFab(lev);

        for (amrex::MFIter mfi(field_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            values.resize(ncells);
            HYPRE_SStructVectorGetBoxValues(x, lev, lo, hi, variable, values.data());

            const auto cval_arr  = cval_mf.const_array(mfi);
            auto       field_arr = field_mf.array(mfi);

            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (cval_arr(ii, jj, kk) >= 0)
                    field_arr(ii, jj, kk) = values[cnt];
                ++cnt;
            }
        }
    }
#else
    HYPRE_SStructVectorGetBoxValues(x, 0, ilower, iupper, variable, values.data());

    count = 0;
    KJILOOP
    {
        if (cval(i, j, k) >= 0)
            f(i, j, k) = values[count];
        ++count;
    }
#endif
}
