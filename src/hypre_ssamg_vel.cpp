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
// Multi-level (nlevs>1): the C-F coupling is handled exactly as for pressure, and for the same
// reason -- the graph carries non-stencil C-F entries (built in make_grid_7p) that MUST be
// filled or HYPRE sees uninitialised couplings. The coupling is entirely in the MATRIX: unlike
// pressure there is no separate gradient/velocity corrector (diffusion is a direct implicit
// solve, not a projection). amr_cf_coefficients is reused verbatim -- it is operator-agnostic,
// rescaling whatever fine-spacing flux idiff2_FS* wrote in the interface slot to the true C-F
// centre distance d_cf=0.5(dx_f+dx_c) (fine scale 2/(1+r)) and making the coarse->fine coupling
// volume-conservative. Covered coarse cells become identity rows; V_lev weighting keeps the
// cross-level operator consistent with that conservative coupling. NO nullspace projection /
// reference pin: the CPOR/(alpha*dt) diagonal makes the Helmholtz operator non-singular.
//
// CAVEAT: the nlevs>1 solve reuses create_solver's PCG+BoomerAMG path, which assumes symmetry.
// The volume-weighted diffusion operator is symmetric on a uniform grid; on a stretched grid the
// interior stencil is only approximately symmetric, so a stretched-grid multi-level run may want
// a GMRES outer solver instead. Single-level (SSAMG/GMRES) is unaffected.
// ====================================================================================

void hypre_ssamg::start_solver123(lexer *p, fdm *a, ghostcell *pgc, field &f, int var)
{
    p->solveriter = 0;

    #if USE_AMREX
    // Rebuild on a regrid, and also whenever the assembled grid's part count no longer
    // matches p->nlevs. p->changed is a one-shot flag that lexer::regrid CLEARS on entry:
    // the init hierarchy build (driver_ini_cfd) runs regrid in a loop whose last, no-growth
    // pass resets it to false, so the level 1->2 creation that happened in the previous pass
    // is never seen here -- the solver objects are constructed by logic_cfd() BEFORE that
    // loop, at nlevs==1. fill_matrix4/fill_matrix_vel then loop lev<p->nlevs and call
    // HYPRE_SStructMatrixSetBoxValues on a part that was never given any extents -> segfault.
    if(p->changed || numparts != p->nlevs)
    {
        destroy_grid();
        make_grid_7p(p, a, pgc);
    }
    #endif

    // Rebuild the solver for every solve. The matrix values change every step, so solve() has
    // to run HYPRE_SStructGMRESSetup/SSAMGSetup each time -- and hypre's SStruct setup allocates
    // a complete new hierarchy on every call (SSAMGSetup -> ComputeRAP / MatvecSetup ->
    // hypre_StructMatrixResize) without releasing the one from the previous call. There is no
    // "unsetup", only Destroy, so a solver object kept across steps leaks one full multigrid
    // hierarchy per solve (measured ~25 MB/step at 240x2x240 with 3 RK3 pressure solves, >5 GB
    // after ~200 steps). Create/Destroy are trivial next to the Setup that has to run anyway.
    bool rebuild_solver = true;

    #if USE_AMREX
    // nlevs>1 takes the ParCSR GMRES + BoomerAMG path, which does NOT leak: hypre_BoomerAMGSetup
    // frees its own previous hierarchy. There the rebuild only has to follow the operator's
    // identity -- grid_rebuilt (make_grid_7p destroyed A/b/x, and a solver kept across that still
    // points at freed matrix data -> segfault in hypre_BoomerAMGCycle) or a change in level count.
    if(p->nlevs > 1)
        rebuild_solver = (!solver_created || created_nlevs != p->nlevs || grid_rebuilt);
    #endif

    if(rebuild_solver)
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
    // Correct the C-F interface coefficients in a->M and fill cf_links[].coeff for the diffusion
    // operator (rescale fine-spacing flux to d_cf, drop the stencil slot, make coarse->fine
    // conservative). Operator-agnostic; returns immediately for nlevs<=1.
    amr_cf_coefficients(p, a, pgc, cval);

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto&       cval_mf  = cval.GetMultiFab(lev);
        const auto& field_mf = f.GetMultiFab(lev);
        const auto& cover_mf = p->amr_cell_mf[lev];

        const auto V_lev = p->amrex_geometry[lev].CellSizeArray()[0]
                         * (p->j_dir ? p->amrex_geometry[lev].CellSizeArray()[1] : 1.0)
                         * p->amrex_geometry[lev].CellSizeArray()[2];

        // Covered coarse cells (refined footprint under a finer patch) are not part of the
        // composite solution -- the fine level is authoritative there. amr_cf_coefficients moved
        // the uncovered coarse cell's coupling into the patch, but the covered neighbour still
        // points back one-sidedly; emit covered cells as identity rows (x=b=0) so they decouple.
        // amr_cell_mf is 1 uncovered / 0 covered; the finest level is a blanket 0.
        const bool has_cover = (lev < p->nlevs - 1);

        for (amrex::MFIter mfi(cval_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto cval_arr  = cval_mf.const_array(mfi);
            const auto field_arr = field_mf.const_array(mfi);
            const auto cover_arr = cover_mf.const_array(mfi);

            auto is_solved = [&](int ii, int jj, int kk)
            {
                return cval_arr(ii, jj, kk) >= 0
                    && !(has_cover && cover_arr(ii, jj, kk) == 0);
            };

            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            // Matrix A: 7 stencil entries per cell (order p, s(x-), n(x+), e(y-), w(y+), b(z-), t(z+)).
            // Interface stencil slots were already zeroed (diagonal corrected) in
            // amr_cf_coefficients; those couplings are emitted below as non-stencil entries.
            values.resize(ncells * 7);
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (is_solved(ii, jj, kk))
                {
                    const int n = cval_arr(ii, jj, kk);
                    values[cnt++] = a->M.p[n] * V_lev;
                    values[cnt++] = a->M.s[n] * V_lev;
                    values[cnt++] = a->M.n[n] * V_lev;
                    values[cnt++] = (ydir ? a->M.e[n] : 0.0) * V_lev;
                    values[cnt++] = (ydir ? a->M.w[n] : 0.0) * V_lev;
                    values[cnt++] = a->M.b[n] * V_lev;
                    values[cnt++] = a->M.t[n] * V_lev;
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
                values[cnt++] = is_solved(ii, jj, kk) ? field_arr(ii, jj, kk) : 0.0;
            HYPRE_SStructVectorSetBoxValues(x, lev, lo, hi, variable, values.data());

            // RHS b (no nullspace projection -- the Helmholtz operator is non-singular).
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
                values[cnt++] = is_solved(ii, jj, kk) ? a->rhsvec.V[cval_arr(ii, jj, kk)] * V_lev : 0.0;
            HYPRE_SStructVectorSetBoxValues(b, lev, lo, hi, variable, values.data());
        }
    }

    // Coarse-fine couplings: replay the non-stencil entries recorded when the graph was built.
    // The .entry slot and per-cell ordering must match make_grid_7p exactly; .coeff was set in
    // amr_cf_coefficients. Only the owner of the "from" cell recorded the link.
    for (const auto& L : cf_links)
    {
        int idx[3] = {L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]};
        int ent    = L.entry;
        const int lev = L.from_part;
        const auto V_lev = p->amrex_geometry[lev].CellSizeArray()[0]
                         * (p->j_dir ? p->amrex_geometry[lev].CellSizeArray()[1] : 1.0)
                         * p->amrex_geometry[lev].CellSizeArray()[2];
        double val = L.coeff * V_lev;
        HYPRE_SStructMatrixSetValues(A, L.from_part, idx, variable, 1, &ent, &val);
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
        const auto& cover_mf = p->amr_cell_mf[lev];
        const bool has_cover = (lev < p->nlevs - 1);

        for (amrex::MFIter mfi(field_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            values.resize(ncells);
            HYPRE_SStructVectorGetBoxValues(x, lev, lo, hi, variable, values.data());

            const auto cval_arr  = cval_mf.const_array(mfi);
            const auto cover_arr = cover_mf.const_array(mfi);
            auto       field_arr = field_mf.array(mfi);

            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                // Only write solved cells; covered coarse cells (identity x=0) keep their
                // pre-solve value (refreshed by average_down), matching fill_matrix_vel's mask.
                if (cval_arr(ii, jj, kk) >= 0
                    && !(has_cover && cover_arr(ii, jj, kk) == 0))
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
