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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"hypre_ssamg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint4.h"

void hypre_ssamg::fill_matrix4(lexer* p, fdm* a, ghostcell* pgc, field& f)
{
    fieldint4 cval4(p);

    count = 0;
    LOOP
    {
        cval4(i, j, k) = count;
        ++count;
    }

    nentries = 7;
    for (j = 0; j < nentries; j++)
        stencil_indices[j] = j;

#if USE_AMREX
    // Build the coarsened bounding box of each fine patch so we can detect
    // coarse cells whose stencil neighbours fall inside the covered region.
    // A coarse stencil entry pointing into the covered region must be zeroed:
    // HYPRE adds the correct inter-part coupling via the AMR machinery and a
    // non-zero stencil entry there would double-count the flux.
    std::vector<amrex::Box> covered_boxes;
    if (p->nlevs > 1)
    {
        amrex::IntVect rf_iv(p->ref_vec[0], p->ref_vec[1], p->ref_vec[2]);
        for (int b = 0; b < (int)p->amrex_box_array[1].size(); ++b)
            covered_boxes.push_back(amrex::coarsen(p->amrex_box_array[1][b], rf_iv));
    }

    auto is_covered = [&](int ci, int cj, int ck) -> bool {
        for (const auto& cb : covered_boxes)
            if (cb.contains(amrex::IntVect(ci, cj, ck))) return true;
        return false;
    };

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto& cval4_mf      = cval4.GetMultiFab(lev);
        const auto& flag4_mf = p->flag4.GetMultiFab(lev);
        const auto& field_mf = f.GetMultiFab(lev);

        for (amrex::MFIter mfi(cval4_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto cval_arr  = cval4_mf.array(mfi);
            const auto flag_arr  = flag4_mf.const_array(mfi);
            const auto field_arr = field_mf.const_array(mfi);

            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            // Matrix A: 7 stencil entries per cell
            values.resize(ncells * 7);
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (flag_arr(ii, jj, kk) >= AIR_FLAG)
                {
                    int n = cval_arr(ii, jj, kk);
                    double vp = a->M.p[n];
                    double vs = a->M.s[n];
                    double vn = a->M.n[n];
                    double ve = a->M.e[n];
                    double vw = a->M.w[n];
                    double vb = a->M.b[n];
                    double vt = a->M.t[n];

                    // On the coarse level (lev==0), zero any off-diagonal entry
                    // whose neighbour falls inside the covered coarse region.
                    // The removed flux is already represented by the inter-part
                    // graph entries added in amr_graph_entries().
                    if (lev == 0)
                    {
                        if (is_covered(ii-1, jj, kk)) { vp -= vs; vs = 0.0; }
                        if (is_covered(ii+1, jj, kk)) { vp -= vn; vn = 0.0; }
                        if (is_covered(ii, jj-1, kk)) { vp -= ve; ve = 0.0; }
                        if (is_covered(ii, jj+1, kk)) { vp -= vw; vw = 0.0; }
                        if (is_covered(ii, jj, kk-1)) { vp -= vb; vb = 0.0; }
                        if (is_covered(ii, jj, kk+1)) { vp -= vt; vt = 0.0; }
                    }

                    values[cnt++] = vp;
                    values[cnt++] = vs;
                    values[cnt++] = vn;
                    values[cnt++] = ve;
                    values[cnt++] = vw;
                    values[cnt++] = vb;
                    values[cnt++] = vt;
                }
                else
                {
                    values[cnt++] = 1.0;
                    values[cnt++] = 0.0;
                    values[cnt++] = 0.0;
                    values[cnt++] = 0.0;
                    values[cnt++] = 0.0;
                    values[cnt++] = 0.0;
                    values[cnt++] = 0.0;
                }
            }
            HYPRE_SStructMatrixSetBoxValues(A, lev, lo, hi, variable, nentries, stencil_indices, values.data());

            // Initial-guess vector x: 1 entry per cell
            values.resize(ncells);
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                values[cnt++] = (flag_arr(ii, jj, kk) >= AIR_FLAG) ? field_arr(ii, jj, kk) : 0.0;
            }
            HYPRE_SStructVectorSetBoxValues(x, lev, lo, hi, variable, values.data());

            // RHS vector b
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (flag_arr(ii, jj, kk) >= AIR_FLAG)
                    values[cnt++] = a->rhsvec.V[cval_arr(ii, jj, kk)];
                else
                    values[cnt++] = 0.0;
            }
            HYPRE_SStructVectorSetBoxValues(b, lev, lo, hi, variable, values.data());
        }
    }

    HYPRE_SStructMatrixAssemble(A);
    HYPRE_SStructVectorAssemble(x);
    HYPRE_SStructVectorAssemble(b);

#else
    values.resize(p->knox * p->knoy * p->knoz * 7);

    // fill matrix coefficients
    count = 0;
    KJILOOP
    {
        PFLUIDCHECK
        {
            n = cval4(i, j, k);

            values[count] = a->M.p[n]; ++count;
            values[count] = a->M.s[n]; ++count;
            values[count] = a->M.n[n]; ++count;
            values[count] = a->M.e[n]; ++count;
            values[count] = a->M.w[n]; ++count;
            values[count] = a->M.b[n]; ++count;
            values[count] = a->M.t[n]; ++count;
        }
        SFLUIDCHECK
        {
            values[count] = 1.0; ++count;
            values[count] = 0.0; ++count;
            values[count] = 0.0; ++count;
            values[count] = 0.0; ++count;
            values[count] = 0.0; ++count;
            values[count] = 0.0; ++count;
            values[count] = 0.0; ++count;
        }
    }

    HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, variable, nentries, stencil_indices, values.data());
    HYPRE_SStructMatrixAssemble(A);

    // initial guess vector x
    count = 0;
    KJILOOP
    {
        PFLUIDCHECK
        values[count] = f(i, j, k);
        SFLUIDCHECK
        values[count] = 0.0;

        ++count;
    }

    HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, variable, values.data());
    HYPRE_SStructVectorAssemble(x);

    // RHS vector b
    count = 0;
    KJILOOP
    {
        PFLUIDCHECK
        {
            n = cval4(i, j, k);
            values[count] = a->rhsvec.V[n];
        }
        SFLUIDCHECK
        values[count] = 0.0;

        ++count;
    }

    HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, variable, values.data());
    HYPRE_SStructVectorAssemble(b);
#endif
}

void hypre_ssamg::fillbackvec4(lexer* p, field& f, int var)
{
#if USE_AMREX
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto& field_mf       = f.GetMultiFab(lev);
        const auto& flag4_mf = p->flag4.GetMultiFab(lev);

        for (amrex::MFIter mfi(field_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            values.resize(ncells);
            HYPRE_SStructVectorGetBoxValues(x, lev, lo, hi, variable, values.data());

            const auto flag_arr  = flag4_mf.const_array(mfi);
            auto       field_arr = field_mf.array(mfi);

            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (flag_arr(ii, jj, kk) >= AIR_FLAG)
                    field_arr(ii, jj, kk) = values[cnt];
                ++cnt;
            }
        }
    }
#else
    HYPRE_SStructVectorGetBoxValues(x, part, ilower, iupper, variable, values.data());

    count = 0;
    KJILOOP
    {
        PFLUIDCHECK
        f(i, j, k) = values[count];

        ++count;
    }
#endif
}
