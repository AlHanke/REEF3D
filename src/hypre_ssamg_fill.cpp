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

#include "hypre_ssamg.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "fieldint4.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>

void hypre_ssamg::fill_matrix4(lexer* p, fdm* a, ghostcell* pgc, field& f)
{
    fieldint4 cval4(p);

    count = 0;
    LOOP
    {
        cval4(i, j, k) = count;
        ++count;
    }

    nentries = stencil_size;
    for (j = 0; j < nentries; j++)
        stencil_indices[j] = j;

#if USE_AMREX
    // Correct the coarse-fine interface coefficients in a->M and fill cf_links[].coeff.
    // This rescales the (fine-spacing) interface entries poisson_pcorr already wrote to
    // the true C-F centre distance, redistributes the coarse flux over the fine
    // sub-faces, and zeros the corresponding stencil slots so the transfer below emits
    // them only as non-stencil entries.
    amr_cf_coefficients(p, a, pgc, cval4);

    // // Rowsum sanity check: for every fluid cell, the sum of its 7-point stencil
    // // entries plus all cf_link coefficients that originate from it must equal the
    // // stencil-only rowsum poisson_pcorr produced (which is zero for interior cells
    // // and equals the wall-BC coefficient for cells adjacent to domain boundaries).
    // // A non-zero residual on interior C-F cells indicates a missing or double-counted
    // // C-F coupling in amr_cf_coefficients.
    // if (p->nlevs > 1)
    // {
    //     // Per-row accumulator for cf_link coefficients (indexed by cval4).
    //     std::vector<double> cfsum(p->veclength, 0.0);
    //     for (int lev = 0; lev < p->nlevs; ++lev)
    //     {
    //         const auto& cmf = cval4.GetMultiFab(lev);
    //         for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
    //         {
    //             const amrex::Box& bx = mfi.validbox();
    //             const auto ca = cmf.const_array(mfi);
    //             for (const auto& L : cf_links)
    //             {
    //                 if (L.from_part != lev) continue;
    //                 const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
    //                 if (!bx.contains(fiv)) continue;
    //                 cfsum[ca(fiv)] += L.coeff;
    //             }
    //         }
    //     }

    //     double maxrs = 0.0;
    //     for (int lev = 0; lev < p->nlevs; ++lev)
    //     {
    //         const auto& cmf  = cval4.GetMultiFab(lev);
    //         const auto& f4mf = p->flag4.GetMultiFab(lev);
    //         for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
    //         {
    //             const amrex::Box& bx = mfi.validbox();
    //             const auto ca = cmf.const_array(mfi);
    //             const auto fa = f4mf.const_array(mfi);
    //             amrex::LoopOnCpu(bx, [&](int ii, int jj, int kk) noexcept
    //             {
    //                 if (fa(ii, jj, kk) < AIR_FLAG) return;
    //                 const int n = ca(ii, jj, kk);
    //                 const double rs = a->M.p[n] + a->M.n[n] + a->M.s[n]
    //                                 + a->M.w[n] + a->M.e[n]
    //                                 + a->M.t[n] + a->M.b[n] + cfsum[n];
    //                 maxrs = std::max(maxrs, std::fabs(rs));
    //             });
    //         }
    //     }
    //     const double g_maxrs = pgc->globalmax(maxrs);
    //     if (p->mpirank == 0)
    //         std::cout << "\n  [rowsum] max|stencil+cflink| = " << g_maxrs
    //                   << "  (interior ~0, wall-BC cells show wall coeff)" << std::endl;
    // }

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
                    // Interface stencil slots were already zeroed (and the diagonal
                    // corrected) in amr_cf_coefficients; those couplings are emitted
                    // separately as non-stencil entries from cf_links.
                    int n = cval_arr(ii, jj, kk);
                    values[cnt++] = a->M.p[n];
                    values[cnt++] = a->M.s[n];
                    values[cnt++] = a->M.n[n];
                    values[cnt++] = a->M.e[n];
                    values[cnt++] = a->M.w[n];
                    values[cnt++] = a->M.b[n];
                    values[cnt++] = a->M.t[n];
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

    // Coarse-fine couplings: replay the non-stencil entries recorded when the graph
    // was built. The .entry slot (stencil_size + k) and per-cell ordering must match
    // make_grid_7p exactly; .coeff is set during matrix preparation. Only the process
    // that owns the "from" cell recorded the link, so no ownership re-check is needed.
    for (const auto& L : cf_links)
    {
        int    idx[3] = {L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]};
        int    ent    = L.entry;
        double val    = L.coeff;
        HYPRE_SStructMatrixSetValues(A, L.from_part, idx, variable, 1, &ent, &val);
    }

    HYPRE_SStructMatrixAssemble(A);
    HYPRE_SStructVectorAssemble(x);
    HYPRE_SStructVectorAssemble(b);

    // Multi-level: hand the assembled ParCSR operator/vectors to the BiCGSTAB+BoomerAMG
    // path. The ParCSR objects share storage with the SStruct ones, so the solution
    // written into par_x is read back through x in fillbackvec4.
    if (p->nlevs > 1)
    {
        HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
        HYPRE_SStructVectorGetObject(b, (void**) &par_b);
        HYPRE_SStructVectorGetObject(x, (void**) &par_x);
    }

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

    HYPRE_SStructMatrixSetBoxValues(A, 0, ilower, iupper, variable, nentries, stencil_indices, values.data());
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

    HYPRE_SStructVectorSetBoxValues(x, 0, ilower, iupper, variable, values.data());
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

    HYPRE_SStructVectorSetBoxValues(b, 0, ilower, iupper, variable, values.data());
    HYPRE_SStructVectorAssemble(b);
#endif

    // Operator pre-checks on the fully assembled A (and, for nlevs>1, the extracted
    // ParCSR object). Off by default; export REEF_MAT_CHECK=1 to enable.
    if (std::getenv("REEF_MAT_CHECK"))
        validate_operator(p, a, pgc);
}

void hypre_ssamg::amr_cf_coefficients(lexer* p, fdm* a, ghostcell* pgc, fieldint4& cval4)
{
#if USE_AMREX
    // Guard on nlevs only (global): the Allgather below is collective, so a rank owning no
    // C-F links must still reach it (its loops simply do nothing).
    if (p->nlevs <= 1) return;

    const auto& rr = p->ref_vec;

    // matrix_diag slot for a face given its axis (0:x,1:y,2:z) and side, matching the
    // 7-point stencil order in make_grid_7p: x-(s)/x+(n), y-(e)/y+(w), z-(b)/z+(t).
    auto slot = [&](int n, int axis, bool high) -> double&
    {
        if (axis == 0) return high ? a->M.n[n] : a->M.s[n];
        if (axis == 1) return high ? a->M.w[n] : a->M.e[n];
        return                high ? a->M.t[n] : a->M.b[n];
    };

    // Pass 1: for every coupling resolve the global row index n (held in cval4, only
    // addressable through the owning level's MFIter array), the interface axis/side in
    // the from-cell's own index space, and the original fine-spacing face coefficient
    // poisson_pcorr wrote there.
    struct rec { int n, axis; bool high, is_fine; double old; };
    std::vector<rec> info(cf_links.size(), {-1, -1, false, false, 0.0});

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        const auto& cval_mf = cval4.GetMultiFab(lev);

        for (amrex::MFIter mfi(cval_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto cval_arr  = cval_mf.const_array(mfi);

            for (int id = 0; id < (int)cf_links.size(); ++id)
            {
                const auto& L = cf_links[id];
                if (L.from_part != lev) continue;

                const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
                if (!bx.contains(fiv)) continue;

                const bool          is_fine = (L.from_part > L.to_part);
                const amrex::IntVect tiv(L.to_ijk[0], L.to_ijk[1], L.to_ijk[2]);

                // Interface normal/side expressed in the from-cell's own index space.
                int  axis = 0;
                bool high = false;
                if (is_fine)
                {
                    const amrex::IntVect cf = amrex::coarsen(fiv, rr);
                    for (int d = 0; d < 3; ++d) if (cf[d] != tiv[d]) axis = d;
                    high = (tiv[axis] > cf[axis]);
                }
                else
                {
                    const amrex::IntVect cf = amrex::coarsen(tiv, rr);
                    for (int d = 0; d < 3; ++d) if (cf[d] != fiv[d]) axis = d;
                    high = (cf[axis] > fiv[axis]);
                }

                const int n = cval_arr(fiv);
                info[id] = {n, axis, high, is_fine, slot(n, axis, high)};
            }
        }
    }

    // Pass 2: couplings sharing a (row, axis, side) belong to one physical C-F face.
    // Per face: rescale poisson's fine-spacing entry to the true C-F centre distance
    // dcf = 0.5*(dx_f+dx_c) -- a pure function of the refinement ratio r:
    //   fine side  dx_f/dcf = 2/(1+r),   coarse side  dx_c/dcf = 2r/(1+r).
    // The diagonal is shifted from the old to the corrected flux once, the stencil slot
    // is dropped, and the corrected flux is split equally over the face's couplings.
    std::map<std::array<int,3>, std::vector<int>> faces; // (n,axis,side) -> link ids
    for (int id = 0; id < (int)cf_links.size(); ++id)
    {
        if (info[id].n < 0) continue; // defensive: not owned here
        faces[{info[id].n, info[id].axis, info[id].high ? 1 : 0}].push_back(id);
    }

    // Pass 2a -- FINE->coarse faces: rescale poisson's fine-spacing entry to the C-F centre
    // distance dcf (fine scale 2/(1+r)). The result embeds the FINE-side face density.
    for (auto& [key, ids] : faces)
    {
        const rec& r0 = info[ids[0]];
        if (!r0.is_fine) continue;
        const double scale = 2.0 / (1.0 + double(rr[r0.axis]));
        const double new_total = r0.old * scale;          // corrected fine flux (<= 0)
        a->M.p[r0.n] += (r0.old - new_total);             // diagonal: T_old -> T_cf
        slot(r0.n, r0.axis, r0.high) = 0.0;               // drop the fine-spacing entry
        const double per = new_total / double(ids.size());
        for (int id : ids) cf_links[id].coeff = per;
    }

    // Gather every fine->coarse coupling globally: the coarse cell of a C-F face is in
    // general on a different MPI rank than its covering fine cells, so it cannot read the
    // fine coefficients locally. Pack [fine_ijk(3), coarse_ijk(3), coeff] and Allgather.
    std::map<std::array<int,6>,double> gfc;
    {
        std::vector<double> loc;
        for (const auto& L : cf_links)
            if (L.from_part > L.to_part) // fine->coarse
                loc.insert(loc.end(), {double(L.from_ijk[0]),double(L.from_ijk[1]),double(L.from_ijk[2]),
                                       double(L.to_ijk[0]),  double(L.to_ijk[1]),  double(L.to_ijk[2]),
                                       L.coeff});
        int ncomm; MPI_Comm_size(pgc->mpi_comm,&ncomm);
        int ncnt=(int)loc.size();
        std::vector<int> counts(ncomm), displs(ncomm);
        MPI_Allgather(&ncnt,1,MPI_INT,counts.data(),1,MPI_INT,pgc->mpi_comm);
        int tot=0; for(int r=0;r<ncomm;++r){ displs[r]=tot; tot+=counts[r]; }
        std::vector<double> all(tot);
        MPI_Allgatherv(loc.data(),ncnt,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,pgc->mpi_comm);
        for (int o=0; o+6<tot; o+=7)
            gfc[{int(all[o]),int(all[o+1]),int(all[o+2]),int(all[o+3]),int(all[o+4]),int(all[o+5])}] = all[o+6];
    }

    // Pass 2b -- COARSE->fine faces: make the coupling conservative. Setting
    // coeff_cf(F) = coeff_fc(F) / (vol_c/vol_f) gives M_cf*vol_c = M_fc*vol_f exactly, so
    // the C-F flux is single-valued (the coarse face density becomes the fine average).
    const double volratio = double(rr[0])*double(rr[1])*double(rr[2]);
    for (auto& [key, ids] : faces)
    {
        const rec& r0 = info[ids[0]];
        if (r0.is_fine) continue;
        double new_total = 0.0;
        for (int id : ids)
        {
            const auto& L = cf_links[id]; // from=coarse, to=fine
            const std::array<int,6> k2{L.to_ijk[0],L.to_ijk[1],L.to_ijk[2],
                                       L.from_ijk[0],L.from_ijk[1],L.from_ijk[2]};
            const auto it = gfc.find(k2);
            const double coeff = (it!=gfc.end() ? it->second : 0.0) / volratio;
            cf_links[id].coeff = coeff;
            new_total += coeff;
        }
        a->M.p[r0.n] += (r0.old - new_total); // diagonal: T_old -> conservative T_cf
        slot(r0.n, r0.axis, r0.high) = 0.0;
    }
#endif
}

void hypre_ssamg::cf_velocity_correction(lexer* p, fdm* a, ghostcell* pgc,
                                         field& pcorr, field& u, field& v, field& w,
                                         double alpha)
{
#if USE_AMREX
    // C-F consistent gradient correction (the G side of the composite projection). The matrix
    // coupled each fine interface cell to the real coarse cell over the C-F centre distance
    // dcf = 0.5*(dx_f+dx_c); the interior ucorr/vcorr/wcorr instead used the interpolated
    // start4 ghost at fine spacing. Per fine->coarse link, replace the interior face correction
    // with the matrix-consistent one so the fine cell's discrete divergence D.(1/rho G pcorr)
    // matches its matrix row L.pcorr at the C-F face.
    //
    // The coefficient is taken from the stored matrix coupling, so density/porosity are
    // identical to those in L (no separate roface evaluation):
    //   |new_total| = |L.coeff|              (the rescaled fine-side C-F entry, over dcf)
    //   |old|       = |new_total|*(1+r)/2    (undo the fine-side scale 2/(1+r) -> over dx_f)
    if (p->nlevs <= 1) return;

    const auto& rr = p->ref_vec;

    // --- TEMPORARY localization probe (REEF_CF_DIV), remove once the C-F ring is fixed ---
    // Measures the post-correction discrete divergence of each fine interface cell, split by
    // link type (high/low/both) and lateral wall distance, and reports the global worst cell.
    const bool dbg = (std::getenv("REEF_CF_DIV") != nullptr);
    double g_hi = 0.0, g_lo = 0.0, g_both = 0.0;   // max|D.u| by link type (interface cells)
    double g_wall = 0.0, g_inner = 0.0;             // max|D.u| by lateral wall distance (<=1 cell)
    double l_worst = 0.0; int l_lev = -1, l_ijk[3] = {-1,-1,-1}; bool l_hi = false, l_lo = false;
    int l_fc = 0, l_fzm = 0, l_fzp = 0, l_fxm = 0, l_fxp = 0; // flag4 at worst cell + z/x neighbours
    double g_wall_noif = 0.0;                        // max|D.u| at lateral-wall NON-interface fluid cells
    double n_worst = 0.0; int n_lev = -1, n_ijk[3] = {-1,-1,-1};
    double g_cwall = 0.0;                            // coarse(lev0) lateral-wall, EXCLUDING covered cells
    double c_worst = 0.0; int c_ijk[3] = {-1,-1,-1};

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        // Coarse pcorr local to each fine box on its coarsened layout, one ghost so the coarse
        // cell across the C-F interface (the matrix's "to" cell) is addressable.
        const auto& crse_mf = pcorr.GetMultiFab(lev - 1);
        auto&       fine_mf = pcorr.GetMultiFab(lev);

        amrex::BoxArray cfba = amrex::coarsen(fine_mf.boxArray(), rr);
        amrex::MultiFab  coarse_on_fine(cfba, fine_mf.DistributionMap(), 1, 1);
        coarse_on_fine.setVal(0.0);
        coarse_on_fine.ParallelCopy(crse_mf, 0, 0, 1,
                                    amrex::IntVect(0), amrex::IntVect(1),
                                    p->amrex_geometry[lev - 1].periodicity());

        auto& u_mf = u.GetMultiFab(lev);
        auto& v_mf = v.GetMultiFab(lev);
        auto& w_mf = w.GetMultiFab(lev);

        const double dx[3] = {p->amrex_geometry[lev].CellSize(0),
                              p->amrex_geometry[lev].CellSize(1),
                              p->amrex_geometry[lev].CellSize(2)};

        const amrex::Box dombx = p->amrex_geometry[lev].Domain();
        const int domlo[3] = {dombx.smallEnd(0), dombx.smallEnd(1), dombx.smallEnd(2)};
        const int domhi[3] = {dombx.bigEnd(0),   dombx.bigEnd(1),   dombx.bigEnd(2)};

        for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();

            const auto pf_arr = fine_mf.const_array(mfi);        // fine pcorr (+ start4 ghost)
            const auto pc_arr = coarse_on_fine.const_array(mfi); // coarse pcorr, coarse indices
            auto       u_arr  = u_mf.array(mfi);
            auto       v_arr  = v_mf.array(mfi);
            auto       w_arr  = w_mf.array(mfi);

            for (const auto& L : cf_links)
            {
                if (L.from_part != lev)         continue; // fine cell lives on this level
                if (L.to_part   >= L.from_part) continue; // fine->coarse couplings only

                const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
                if (!bx.contains(fiv)) continue;

                const amrex::IntVect tiv(L.to_ijk[0], L.to_ijk[1], L.to_ijk[2]);
                const amrex::IntVect cf = amrex::coarsen(fiv, rr);

                // Interface axis/side in the fine cell's own index space.
                int axis = 0;
                for (int d = 0; d < 3; ++d) if (cf[d] != tiv[d]) axis = d;
                const bool high = (tiv[axis] > cf[axis]);

                const double absnew = std::fabs(L.coeff);                      // over dcf
                const double absold = absnew * 0.5 * (1.0 + double(rr[axis]));  // over dx_f

                const double pf = pf_arr(fiv);
                const double pc = pc_arr(tiv);

                amrex::IntVect e(0, 0, 0); e[axis] = 1;

                if (high)
                {
                    // Owned high face: undo the interior correction (interpolated ghost at fine
                    // spacing) and apply the C-F one (real coarse pcorr over dcf).
                    const double pg = pf_arr(fiv + e); // interpolated start4 ghost
                    const double delta = alpha * p->dt * dx[axis]
                                       * (absnew * (pc - pf) - absold * (pg - pf));
                    if      (axis == 0) u_arr(fiv) -= delta;
                    else if (axis == 1) v_arr(fiv) -= delta;
                    else                w_arr(fiv) -= delta;
                }
                else
                {
                    // Low face (a ghost the interior loop never touched): full C-F correction.
                    const amrex::IntVect lf = fiv - e;
                    const double delta = alpha * p->dt * dx[axis] * (absnew * (pf - pc));
                    if      (axis == 0) u_arr(lf) -= delta;
                    else if (axis == 1) v_arr(lf) -= delta;
                    else                w_arr(lf) -= delta;
                }
            }

            if (dbg)
            {
                const auto fa = p->flag4.GetMultiFab(lev).const_array(mfi);

                // Unique fine interface cells in this box, with high/low link presence.
                std::map<std::array<int,3>, std::pair<bool,bool>> seen;
                for (const auto& L : cf_links)
                {
                    if (L.from_part != lev || L.to_part >= L.from_part) continue;
                    const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
                    if (!bx.contains(fiv)) continue;
                    const amrex::IntVect cf  = amrex::coarsen(fiv, rr);
                    const amrex::IntVect tiv(L.to_ijk[0], L.to_ijk[1], L.to_ijk[2]);
                    int ax = 0; for (int d = 0; d < 3; ++d) if (cf[d] != tiv[d]) ax = d;
                    auto& hl = seen[{fiv[0], fiv[1], fiv[2]}];
                    if (tiv[ax] > cf[ax]) hl.first = true; else hl.second = true;
                }
                for (auto& [ijk, hl] : seen)
                {
                    const int ci = ijk[0], cj = ijk[1], ck = ijk[2];
                    double d = (u_arr(ci,cj,ck) - u_arr(ci-1,cj,ck)) / dx[0]
                             + (w_arr(ci,cj,ck) - w_arr(ci,cj,ck-1)) / dx[2];
                    if (p->j_dir) d += (v_arr(ci,cj,ck) - v_arr(ci,cj-1,ck)) / dx[1];
                    const double ad = std::fabs(d);

                    if      (hl.first && hl.second) { if (ad > g_both) g_both = ad; }
                    else if (hl.first)              { if (ad > g_hi)   g_hi   = ad; }
                    else                            { if (ad > g_lo)   g_lo   = ad; }

                    int dist = (ci - domlo[0] < domhi[0] - ci) ? ci - domlo[0] : domhi[0] - ci;
                    if (p->j_dir)
                    {
                        const int dy = (cj - domlo[1] < domhi[1] - cj) ? cj - domlo[1] : domhi[1] - cj;
                        if (dy < dist) dist = dy;
                    }
                    if (dist <= 1) { if (ad > g_wall) g_wall = ad; }
                    else           { if (ad > g_inner) g_inner = ad; }

                    if (ad > l_worst)
                    {
                        l_worst = ad; l_lev = lev; l_ijk[0]=ci; l_ijk[1]=cj; l_ijk[2]=ck;
                        l_hi=hl.first; l_lo=hl.second;
                        l_fc  = fa(ci,cj,ck);
                        l_fzm = fa(ci,cj,ck-1); l_fzp = fa(ci,cj,ck+1);
                        l_fxm = fa(ci-1,cj,ck); l_fxp = fa(ci+1,cj,ck);
                    }
                }

                // Disambiguation: divergence at lateral-wall fluid cells that are NOT interface
                // cells. If this matches the interface-wall value, the residual is the general
                // wall-BC treatment; if it is ~0, the error is specific to C-F + wall corners.
                for (int ck = bx.smallEnd(2); ck <= bx.bigEnd(2); ++ck)
                for (int cj = bx.smallEnd(1); cj <= bx.bigEnd(1); ++cj)
                for (int ci = bx.smallEnd(0); ci <= bx.bigEnd(0); ++ci)
                {
                    const bool onwall = (ci==domlo[0] || ci==domhi[0]
                                       || (p->j_dir && (cj==domlo[1] || cj==domhi[1])));
                    if (!onwall) continue;
                    if (fa(ci,cj,ck) < AIR_FLAG) continue;          // fluid cells only
                    if (seen.count({ci,cj,ck})) continue;           // skip interface cells
                    double d = (u_arr(ci,cj,ck) - u_arr(ci-1,cj,ck)) / dx[0]
                             + (w_arr(ci,cj,ck) - w_arr(ci,cj,ck-1)) / dx[2];
                    if (p->j_dir) d += (v_arr(ci,cj,ck) - v_arr(ci,cj-1,ck)) / dx[1];
                    const double ad = std::fabs(d);
                    if (ad > g_wall_noif) g_wall_noif = ad;
                    if (ad > n_worst) { n_worst=ad; n_lev=lev; n_ijk[0]=ci; n_ijk[1]=cj; n_ijk[2]=ck; }
                }
            }
        }
    }

    if (dbg)
    {
        // Coarse (level 0) lateral-wall divergence, EXCLUDING cells under the fine patch, so this
        // samples pure single-level-like coarse wall cells. Nonzero & ~pcorr-scaled here proves the
        // wall-BC error is pre-existing (independent of AMR), just ~4x smaller than fine (T~1/dx^2).
        {
            const double dx0[3] = {p->amrex_geometry[0].CellSize(0),
                                   p->amrex_geometry[0].CellSize(1),
                                   p->amrex_geometry[0].CellSize(2)};
            const amrex::Box dom0 = p->amrex_geometry[0].Domain();
            const int dl0[3] = {dom0.smallEnd(0), dom0.smallEnd(1), dom0.smallEnd(2)};
            const int dh0[3] = {dom0.bigEnd(0),   dom0.bigEnd(1),   dom0.bigEnd(2)};
            const amrex::BoxArray covered = amrex::coarsen(p->amrex_box_array[1], rr);
            auto& u0 = u.GetMultiFab(0);
            for (amrex::MFIter mfi(u0); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const auto ua = u0.const_array(mfi);
                const auto va = v.GetMultiFab(0).const_array(mfi);
                const auto wa = w.GetMultiFab(0).const_array(mfi);
                const auto fa = p->flag4.GetMultiFab(0).const_array(mfi);
                for (int ck = bx.smallEnd(2); ck <= bx.bigEnd(2); ++ck)
                for (int cj = bx.smallEnd(1); cj <= bx.bigEnd(1); ++cj)
                for (int ci = bx.smallEnd(0); ci <= bx.bigEnd(0); ++ci)
                {
                    const bool onwall = (ci==dl0[0] || ci==dh0[0]
                                       || (p->j_dir && (cj==dl0[1] || cj==dh0[1])));
                    if (!onwall) continue;
                    if (fa(ci,cj,ck) < AIR_FLAG) continue;
                    if (covered.contains(amrex::IntVect(ci,cj,ck))) continue;   // skip cells under fine patch
                    double d = (ua(ci,cj,ck) - ua(ci-1,cj,ck)) / dx0[0]
                             + (wa(ci,cj,ck) - wa(ci,cj,ck-1)) / dx0[2];
                    if (p->j_dir) d += (va(ci,cj,ck) - va(ci,cj-1,ck)) / dx0[1];
                    const double ad = std::fabs(d);
                    if (ad > g_cwall)  g_cwall = ad;
                    if (ad > c_worst) { c_worst=ad; c_ijk[0]=ci; c_ijk[1]=cj; c_ijk[2]=ck; }
                }
            }
        }

        const double r_hi    = pgc->globalmax(g_hi);
        const double r_lo    = pgc->globalmax(g_lo);
        const double r_both  = pgc->globalmax(g_both);
        const double r_wall  = pgc->globalmax(g_wall);
        const double r_inner = pgc->globalmax(g_inner);
        const double r_worst = pgc->globalmax(l_worst);
        const double r_wall_noif = pgc->globalmax(g_wall_noif);
        const double r_nworst    = pgc->globalmax(n_worst);
        const double r_cwall     = pgc->globalmax(g_cwall);
        const double r_cworst    = pgc->globalmax(c_worst);
        if (p->mpirank == 0)
        {
            std::cout << "\n  [cf-div] max|D.u| HIGH=" << r_hi << " LOW=" << r_lo
                      << " BOTH=" << r_both << "  | wall(<=1)=" << r_wall
                      << " inner=" << r_inner << std::endl;
            std::cout << "  [cf-div] fine lateral-wall NON-interface max|D.u|=" << r_wall_noif
                      << "  | coarse(lev0) uncovered wall max|D.u|=" << r_cwall << std::endl;
        }
        if (l_lev >= 0 && l_worst == r_worst)   // the rank holding the global worst prints it
            std::cout << "  [cf-div] worst (interface) |D.u|=" << r_worst << " at lev=" << l_lev
                      << " (" << l_ijk[0] << "," << l_ijk[1] << "," << l_ijk[2]
                      << ") hi=" << l_hi << " lo=" << l_lo
                      << "  flag4[c]=" << l_fc
                      << " z-=" << l_fzm << " z+=" << l_fzp
                      << " x-=" << l_fxm << " x+=" << l_fxp
                      << "  (>-1 = water)" << std::endl;
        if (n_lev >= 0 && n_worst == r_nworst)
            std::cout << "  [cf-div] worst (non-if wall) |D.u|=" << r_nworst << " at lev=" << n_lev
                      << " (" << n_ijk[0] << "," << n_ijk[1] << "," << n_ijk[2] << ")" << std::endl;
        if (c_ijk[0] >= 0 && c_worst == r_cworst)
            std::cout << "  [cf-div] worst coarse wall |D.u|=" << r_cworst
                      << " at lev=0 (" << c_ijk[0] << "," << c_ijk[1] << "," << c_ijk[2] << ")" << std::endl;
    }
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
    HYPRE_SStructVectorGetBoxValues(x, 0, ilower, iupper, variable, values.data());

    count = 0;
    KJILOOP
    {
        PFLUIDCHECK
        f(i, j, k) = values[count];

        ++count;
    }
#endif
}
