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
#include "density.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#if USE_AMREX
#include <AMReX_MultiFabUtil.H>
#endif

void hypre_ssamg::fill_matrix4(lexer* p, fdm* a, ghostcell* pgc, field& f)
{
    // Row numbering is owned by poisson_pcorr::start (single source of truth): it stamps
    // a->Mrow(i,j,k) with the LOOP index it used to fill a->M / a->rhsvec. Alias it here so
    // the existing cval4-based indexing (and amr_cf_coefficients) consumes that persistent
    // field directly. Rebuilding a local LOOP numbering here risked drifting from poisson's
    // enumeration whenever the two loop macros disagreed (interior solids, inflow/outflow).
    fieldint4& cval4 = a->Mrow;

    nentries = stencil_size;
    for (j = 0; j < nentries; j++)
        stencil_indices[j] = j;

#if USE_AMREX
    // REEF_CONTAINMENT_CHECK: verify the density transition (|phi| < psi) is fully inside the fine
    // level. The band tags |phi| < 2*psi (ini_psi: psi = F45/3*(dx+dy+dz); band = 2*that), so every
    // UNCOVERED coarse fluid cell should have |phi| >= 2*psi -- and certainly >= psi -- i.e. sit in a
    // pure phase. If the smallest |phi| over uncovered coarse cells drops below psi, the smeared
    // density jump reaches the coarse level (NOT contained): gridding/proper-nesting eroded the
    // patch, the regrid lagged the interface, or the band is too narrow.
    if (std::getenv("REEF_CONTAINMENT_CHECK") && p->nlevs > 1)
    {
        const double psi = p->psi;
        for (int lev = 0; lev < p->nlevs - 1; ++lev)
        {
            const auto& cover = p->amr_cell_mf[lev];
            const auto& phimf = a->phi.GetMultiFab(lev);
            const auto& fmf   = p->flag4.GetMultiFab(lev);
            double loc_min = 1e300; int mi=-1, mj=-1, mk=-1;
            for (amrex::MFIter mfi(cover); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const auto cov = cover.const_array(mfi);
                const auto ph  = phimf.const_array(mfi);
                const auto fa  = fmf.const_array(mfi);
                for (int kk=bx.smallEnd(2); kk<=bx.bigEnd(2); ++kk)
                for (int jj=bx.smallEnd(1); jj<=bx.bigEnd(1); ++jj)
                for (int ii=bx.smallEnd(0); ii<=bx.bigEnd(0); ++ii)
                {
                    if (fa(ii,jj,kk) < AIR_FLAG) continue;   // not fluid
                    if (cov(ii,jj,kk) == 0)      continue;   // covered -> fine authoritative
                    const double aph = std::fabs(ph(ii,jj,kk));
                    if (aph < loc_min) { loc_min=aph; mi=ii; mj=jj; mk=kk; }
                }
            }
            const double g_min = -pgc->globalmax(-loc_min);
            if (p->mpirank == 0)
                std::cout << "  [contain] lev=" << lev << "  psi=" << psi << "  band=" << 2*psi
                          << "  min|phi| over uncovered fluid=" << g_min
                          << (g_min >= psi ? "  -> jump CONTAINED" : "  -> NOT CONTAINED") << std::endl;
            if (loc_min == g_min && mi >= 0)
                std::cout << "    [contain] closest uncovered cell lev=" << lev << " ("
                          << mi << "," << mj << "," << mk << ") |phi|=" << g_min << std::endl;
        }
    }

    // Correct the coarse-fine interface coefficients in a->M and fill cf_links[].coeff.
    // This rescales the (fine-spacing) interface entries poisson_pcorr already wrote to
    // the true C-F centre distance, redistributes the coarse flux over the fine
    // sub-faces, and zeros the corresponding stencil slots so the transfer below emits
    // them only as non-stencil entries.
    amr_cf_coefficients(p, a, pgc, cval4);

    // Rowsum sanity check: for every fluid cell, the sum of its 7-point stencil
    // entries plus all cf_link coefficients that originate from it must equal the
    // stencil-only rowsum poisson_pcorr produced (which is zero for interior cells
    // and equals the wall-BC coefficient for cells adjacent to domain boundaries).
    // A non-zero residual on interior C-F cells indicates a missing or double-counted
    // C-F coupling in amr_cf_coefficients.
    if (p->nlevs > 1)
    {
        // Per-row accumulator for cf_link coefficients (indexed by cval4).
        std::vector<double> cfsum(p->veclength, 0.0);
        for (int lev = 0; lev < p->nlevs; ++lev)
        {
            const auto& cmf = cval4.GetMultiFab(lev);
            for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const auto ca = cmf.const_array(mfi);
                for (const auto& L : cf_links)
                {
                    if (L.from_part != lev) continue;
                    const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
                    if (!bx.contains(fiv)) continue;
                    cfsum[ca(fiv)] += L.coeff;
                }
            }
        }

        // REEF_ROWSUM_PROBE: report the worst |stencil+cflink| rowsum and, at that cell,
        // dump every cf_link originating from it (to_ijk / axis / high / coeff). A duplicated
        // fine->coarse link shows up as >1 link for the same (axis,side) and a rowsum imbalance
        // of ~one face coefficient. Compare single-box vs multi-box to localise the multi-box
        // C-F operator inconsistency.
        if (std::getenv("REEF_ROWSUM_PROBE"))
        {
            // cf_links census (this rank): count fine->coarse links per axis and the max
            // from-index each axis reaches, to reveal whether the x-pos edge (x=79) has links.
            int fc_axis[3]={0,0,0}, fc_maxfrom[3]={-1,-1,-1}, fc_x79=0;
            for (const auto& L : cf_links)
            {
                if (L.to_part >= L.from_part) continue; // fine->coarse only
                if (L.axis>=0 && L.axis<3)
                {
                    ++fc_axis[L.axis];
                    fc_maxfrom[L.axis] = std::max(fc_maxfrom[L.axis], L.from_ijk[L.axis]);
                }
                if (L.axis==0 && L.from_ijk[0]==79) ++fc_x79;
            }
            std::cout << "  [cflinks] rank=" << p->mpirank
                      << "  fine->coarse per-axis x/y/z=" << fc_axis[0] << "/" << fc_axis[1] << "/" << fc_axis[2]
                      << "  max from-index x/y/z=" << fc_maxfrom[0] << "/" << fc_maxfrom[1] << "/" << fc_maxfrom[2]
                      << "  links@x=79: " << fc_x79 << std::endl;

            double maxrs = 0.0; int w_lev=-1, w_ijk[3]={-1,-1,-1};
            double w_stncl=0.0, w_cf=0.0;
            for (int lev = 0; lev < p->nlevs; ++lev)
            {
                const auto& cmf  = cval4.GetMultiFab(lev);
                const auto& f4mf = p->flag4.GetMultiFab(lev);
                for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    const auto ca = cmf.const_array(mfi);
                    const auto fa = f4mf.const_array(mfi);
                    for (int kk=bx.smallEnd(2); kk<=bx.bigEnd(2); ++kk)
                    for (int jj=bx.smallEnd(1); jj<=bx.bigEnd(1); ++jj)
                    for (int ii=bx.smallEnd(0); ii<=bx.bigEnd(0); ++ii)
                    {
                        if (fa(ii,jj,kk) < AIR_FLAG) continue;
                        const int n = ca(ii,jj,kk);
                        const double stncl = a->M.p[n] + a->M.n[n] + a->M.s[n]
                                           + a->M.w[n] + a->M.e[n] + a->M.t[n] + a->M.b[n];
                        const double rs = stncl + cfsum[n];
                        if (std::fabs(rs) > maxrs)
                        {
                            maxrs = std::fabs(rs); w_lev=lev;
                            w_ijk[0]=ii; w_ijk[1]=jj; w_ijk[2]=kk;
                            w_stncl=stncl; w_cf=cfsum[n];
                        }
                    }
                }
            }
            const double g_maxrs = pgc->globalmax(maxrs);
            if (p->mpirank == 0)
                std::cout << "\n  [rowsum] max|stencil+cflink| = " << g_maxrs
                          << "  (interior ~0, wall-BC cells show wall coeff)" << std::endl;
            // The rank holding the global worst dumps its cell + originating cf_links.
            if (w_lev >= 0 && maxrs == g_maxrs)
            {
                int nlinks = 0;
                std::cout << "  [rowsum] worst lev=" << w_lev
                          << " (" << w_ijk[0] << "," << w_ijk[1] << "," << w_ijk[2] << ")"
                          << "  stencilsum=" << w_stncl << "  cflinksum=" << w_cf
                          << "  rank=" << p->mpirank << std::endl;
                for (const auto& L : cf_links)
                {
                    if (L.from_part != w_lev) continue;
                    if (L.from_ijk[0]!=w_ijk[0] || L.from_ijk[1]!=w_ijk[1] || L.from_ijk[2]!=w_ijk[2]) continue;
                    ++nlinks;
                    std::cout << "  [rowsum]   link#" << nlinks
                              << " -> part=" << L.to_part
                              << " to=(" << L.to_ijk[0] << "," << L.to_ijk[1] << "," << L.to_ijk[2] << ")"
                              << " axis=" << L.axis << " high=" << L.high
                              << " coeff=" << L.coeff << " entry=" << L.entry << std::endl;
                }
                std::cout << "  [rowsum]   total cf_links from this cell = " << nlinks << std::endl;
            }
        }
    }

    // RHS nullspace projection (multi-level all-Neumann composite operator). The symmetric
    // operator's nullspace is the DOF-space constant, so a solution exists only if 1.b = 0.
    // We project b -= (1.b)/N over the solved DOFs (1.b is the plain sum of the emitted
    // V-weighted b, since the nullvector is the constant in DOF space). The RHS is already
    // conservative to roundoff, so this is essentially a no-op per solve, but it stops the
    // constant carried in the previous-step initial guess from drifting over a long run and
    // keeps PCG's residual meaningful. Diagnostic printed under REEF_RHS_CHECK.
    const bool rhs_check = (std::getenv("REEF_RHS_CHECK") != nullptr);
    const bool project   = (p->nlevs > 1);
    double proj_mean = 0.0;
    if (project || rhs_check)
    {
        double s = 0.0, sabs = 0.0;
        long   nn = 0;
        for (int lev = 0; lev < p->nlevs; ++lev)
        {
            const auto& cmf = cval4.GetMultiFab(lev);
            const auto& fmf = p->flag4.GetMultiFab(lev);
            const auto& vmf = p->amr_cell_mf[lev];
            const auto Vl = p->amrex_geometry[lev].CellSizeArray()[0] * (p->j_dir ? p->amrex_geometry[lev].CellSizeArray()[1] : 1.0) * p->amrex_geometry[lev].CellSizeArray()[2];
            const bool hc = (lev < p->nlevs - 1);
            for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const auto ca = cmf.const_array(mfi);
                const auto fa = fmf.const_array(mfi);
                const auto va = vmf.const_array(mfi);
                amrex::LoopOnCpu(bx, [&](int ii, int jj, int kk) noexcept
                {
                    if (fa(ii, jj, kk) < AIR_FLAG) return;   // identity row
                    if (hc && va(ii, jj, kk) == 0) return;   // covered identity row
                    const double bval = a->rhsvec.V[ca(ii, jj, kk)] * Vl;
                    s += bval; sabs += std::fabs(bval); ++nn;
                });
            }
        }
        const double g_sum = pgc->globalsum(s);
        const double g_abs = pgc->globalsum(sabs);
        const double g_n   = pgc->globalsum(double(nn));
        proj_mean = (project && g_n > 0.0) ? g_sum / g_n : 0.0;
        if (rhs_check && p->mpirank == 0)
            std::cout << "\n  [rhscheck] 1.b = " << g_sum
                      << "  (rel " << (g_abs > 0.0 ? std::fabs(g_sum) / g_abs : 0.0)
                      << ",  mean/DOF " << (g_n > 0.0 ? g_sum / g_n : 0.0)
                      << ",  N=" << long(g_n) << ")"
                      << (project ? "  -- projected out" : "") << std::endl;
    }

    // REEF_BAND_PIN (diagnostic gate): pin one DOF on the finest level to a homogeneous Dirichlet
    // datum (p=0). The thin adaptive interface band is a near-singular Neumann operator -- the
    // "band-drift" mode (band pressure shifting as a whole relative to the anchored coarse) has a
    // tiny eigenvalue, so GMRES converges the residual but the solution grows ~1e10 along it. One
    // pin on the (connected) band removes that mode. Deterministic pick: the globally-lowest-index
    // fluid cell on the finest level, so the choice is rank- and layout-independent. This confirms
    // the diagnosis; the permanent fix is the covered-cell restriction coupling (no pin). NOTE:
    // assumes the fine region is a single connected component -- a fragmented band needs one pin
    // per component.
    const bool band_pin = (std::getenv("REEF_BAND_PIN") != nullptr) && (p->nlevs > 1);
    int pin_lev = -1, pin_ijk[3] = {-1, -1, -1};
    if (band_pin)
    {
        const int flev = p->nlevs - 1;                 // finest level carries the band
        const auto& cmf = cval4.GetMultiFab(flev);
        const auto& fmf = p->flag4.GetMultiFab(flev);
        double local_key = 1e300; int loc[3] = {-1, -1, -1};
        for (amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto fa = fmf.const_array(mfi);
            for (int kk = bx.smallEnd(2); kk <= bx.bigEnd(2); ++kk)
            for (int jj = bx.smallEnd(1); jj <= bx.bigEnd(1); ++jj)
            for (int ii = bx.smallEnd(0); ii <= bx.bigEnd(0); ++ii)
            {
                if (fa(ii, jj, kk) < AIR_FLAG) continue;
                const double key = (double(ii) * 8192.0 + double(jj)) * 8192.0 + double(kk);
                if (key < local_key) { local_key = key; loc[0]=ii; loc[1]=jj; loc[2]=kk; }
            }
        }
        const double gkey = -pgc->globalmax(-local_key);   // global min key (collective on all ranks)
        if (loc[0] >= 0 && local_key == gkey)
        {
            pin_lev = flev; pin_ijk[0]=loc[0]; pin_ijk[1]=loc[1]; pin_ijk[2]=loc[2];
            std::cout << "  [bandpin] pin lev=" << flev << " (" << loc[0] << "," << loc[1]
                      << "," << loc[2] << ")=0  rank=" << p->mpirank << std::endl;
        }
    }

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        auto& cval4_mf      = cval4.GetMultiFab(lev);
        const auto& flag4_mf = p->flag4.GetMultiFab(lev);
        const auto& field_mf = f.GetMultiFab(lev);
        const auto& cover_mf = p->amr_cell_mf[lev];

        const auto V_lev = p->amrex_geometry[lev].CellSizeArray()[0] * (p->j_dir ? p->amrex_geometry[lev].CellSizeArray()[1] : 1.0) * p->amrex_geometry[lev].CellSizeArray()[2];

        // Covered coarse cells (refined footprint under a finer patch) are NOT part of the
        // composite solution: the fine level is authoritative there. amr_cf_coefficients
        // already dropped the uncovered coarse cell's stencil entry into the patch, but the
        // covered neighbour still points back at it -- a one-sided coupling that breaks
        // symmetry. Emit covered cells as identity rows (x=b=0) so they decouple cleanly.
        // amr_cell_mf is 1 on uncovered / 0 on covered (makeFineMask 1,0); the finest level
        // is a blanket 0, so the mask only marks covered for lev < nlevs-1.
        const bool has_cover = (lev < p->nlevs - 1);

        for (amrex::MFIter mfi(cval4_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto cval_arr  = cval4_mf.array(mfi);
            const auto flag_arr  = flag4_mf.const_array(mfi);
            const auto field_arr = field_mf.const_array(mfi);
            const auto cover_arr = cover_mf.const_array(mfi);

            auto is_solved = [&](int ii, int jj, int kk)
            {
                // Pinned cell -> not solved: emitted as an identity row (x=b=0 => p=0), a single
                // Dirichlet datum that removes the band-drift near-null mode (REEF_BAND_PIN).
                if (band_pin && lev == pin_lev
                    && ii == pin_ijk[0] && jj == pin_ijk[1] && kk == pin_ijk[2]) return false;
                return flag_arr(ii, jj, kk) >= AIR_FLAG
                    && !(has_cover && cover_arr(ii, jj, kk) == 0);
            };

            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const int ncells = bx.numPts();

            // Non-finite pinpoint: report the first solved cell whose matrix row,
            // initial guess, or RHS is Inf/NaN, with level + index + which quantity,
            // so a HYPRE "INFs/NaNs in input" abort can be traced to its source field.
            if (std::getenv("REEF_HYPRE_FINITE_CHECK"))
            {
                for (int kk = lo[2]; kk <= hi[2]; ++kk)
                for (int jj = lo[1]; jj <= hi[1]; ++jj)
                for (int ii = lo[0]; ii <= hi[0]; ++ii)
                {
                    if (!is_solved(ii, jj, kk)) continue;
                    const int n = cval_arr(ii, jj, kk);
                    const double vals[9] = { a->M.p[n], a->M.s[n], a->M.n[n],
                                             a->M.e[n], a->M.w[n], a->M.b[n], a->M.t[n],
                                             field_arr(ii, jj, kk), a->rhsvec.V[n] };
                    const char* nm[9] = {"M.p","M.s","M.n","M.e","M.w","M.b","M.t",
                                         "x(guess)","rhs"};
                    for (int q = 0; q < 9; ++q)
                        if (!std::isfinite(vals[q]))
                        {
                            std::cout << "  [finitecheck] rank " << p->mpirank
                                      << " lev " << lev << " (" << ii << "," << jj << "," << kk
                                      << ") n=" << n << "  " << nm[q] << "=" << vals[q]
                                      << "  flag4=" << flag_arr(ii, jj, kk) << std::endl;
                            break;
                        }
                }
            }

            // Matrix A: 7 stencil entries per cell
            values.resize(ncells * 7);
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (is_solved(ii, jj, kk))
                {
                    // Interface stencil slots were already zeroed (and the diagonal
                    // corrected) in amr_cf_coefficients; those couplings are emitted
                    // separately as non-stencil entries from cf_links.
                    int n = cval_arr(ii, jj, kk);
                    values[cnt++] = a->M.p[n] * V_lev;
                    values[cnt++] = a->M.s[n] * V_lev;
                    values[cnt++] = a->M.n[n] * V_lev;
                    values[cnt++] = a->M.e[n] * V_lev;
                    values[cnt++] = a->M.w[n] * V_lev;
                    values[cnt++] = a->M.b[n] * V_lev;
                    values[cnt++] = a->M.t[n] * V_lev;
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
                values[cnt++] = is_solved(ii, jj, kk) ? field_arr(ii, jj, kk) : 0.0;
            }
            HYPRE_SStructVectorSetBoxValues(x, lev, lo, hi, variable, values.data());

            // RHS vector b
            cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (is_solved(ii, jj, kk))
                    values[cnt++] = a->rhsvec.V[cval_arr(ii, jj, kk)] * V_lev - proj_mean;
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
        // The pinned cell is a pure identity row: zero its outgoing C-F entries too, else the
        // non-stencil couplings would re-populate the row that is_solved emitted as identity.
        const bool pinned = (band_pin && L.from_part == pin_lev
            && L.from_ijk[0]==pin_ijk[0] && L.from_ijk[1]==pin_ijk[1] && L.from_ijk[2]==pin_ijk[2]);

        int    idx[3] = {L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]};
        int    ent    = L.entry;
        const int lev   = L.from_part;
        const auto V_lev = p->amrex_geometry[lev].CellSizeArray()[0] * (p->j_dir ? p->amrex_geometry[lev].CellSizeArray()[1] : 1.0) * p->amrex_geometry[lev].CellSizeArray()[2];
        double val    = pinned ? 0.0 : L.coeff * V_lev;
        HYPRE_SStructMatrixSetValues(A, L.from_part, idx, variable, 1, &ent, &val);
    }

    HYPRE_SStructMatrixAssemble(A);
    HYPRE_SStructVectorAssemble(x);
    HYPRE_SStructVectorAssemble(b);

    // Multi-level: hand the assembled ParCSR operator/vectors to the PCG+BoomerAMG
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
    // ==================================================================================
    // REEF_CF_PROJECTION_GROUP  (grep this tag to find every coupled member)
    // ----------------------------------------------------------------------------------
    // The coarse-fine pressure projection is consistent (assembled operator A == G_v^T W G_v,
    // so the pressure gauge is invariant and NO reference pin / hydrostatic init is required)
    // only if EVERY method in this group encodes the SAME C-F face discretization:
    //   * normal centre distance  d_cf = 0.5*(dx_f + dx_c)   [fine scale 2/(1+r), coarse 2r/(1+r)]
    //   * fine-side face density and porosity. The density now comes from a->rofx/rofy/rofz
    //     (density::update_faces), NOT from a per-call-site pd->roface(): with AMR levels that
    //     array is C-F-synced via amrex::average_down_faces, so the coarse and fine sides of a
    //     C-F face carry one number even when the smoothed Heaviside band straddles it. This
    //     mirrors MLABecLaplacian::averageDownCoeffsToCoarseAmrLevel, which is why amrex_solver
    //     never had the two-valued-density problem.
    //   * conservative volume weighting  (M_cf*V_c == M_fc*V_f, i.e. coeff_cf = coeff_fc/volratio)
    // Change any one of these conventions and you MUST change the others in lockstep, or the
    // fine block floats off its gauge -> spurious velocity / wrong pressure level.
    // See AMR_HYDROSTATIC_FIX_RECORD.md for the derivation and the gauge-check evidence.
    //
    // Members:
    //   1. hypre_ssamg::amr_cf_coefficients     (this file)        matrix C-F coupling            -> defines G_A
    //   2. hypre_ssamg::cf_velocity_correction  (this file)        sole-writer fine C-F vel faces -> G_v (corrector)
    //   3. field_amrex::FillCoarseFineCellGhost (field_amrex.cpp)  C-F ghost ring (predictor/pcorr) -> G_v
    //   4. field_amrex::FillDomainBoundaryImpl  (field_amrex.h)    gcv==41 gate that invokes (3)
    //   5. pjm_corr::start  (gcval_press select) (pjm_corr.cpp)    picks gcv 41 when nlevs>1
    //   6. pjm_corr u/v/wcorr + u/v/wpgrad        (pjm_corr.cpp)    consumers: fine-spacing DXP & a->rof*
    //   7. density::update_faces                  (density.cpp)      sole writer of a->rof*; the C-F sync
    // ==================================================================================

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
                const int  axis = L.axis;
                const bool high = L.high;

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
        if (info[id].n < 0) continue; // defensive: not resolvable in this rank's cval4
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
    // REEF_CF_COARSE_PROBE: does every coarse->fine sub-face find its matching fine->coarse
    // coupling in the Allgathered gfc map? A miss forces coeff=0 -> the coarse cell decouples
    // from the fine region in the matrix (Lp~0) while the reflux still moves the face (D.u!=0),
    // which is exactly the coarse patch-corner residual the adjoint probe localised.
    const bool cprobe = (std::getenv("REEF_CF_COARSE_PROBE") != nullptr);
    long g_hit = 0, g_miss = 0;

    const double volratio = double(rr[0])*double(rr[1])*double(rr[2]);
    for (auto& [key, ids] : faces)
    {
        const rec& r0 = info[ids[0]];
        if (r0.is_fine) continue;
        double new_total = 0.0;
        int nmiss = 0;
        for (int id : ids)
        {
            const auto& L = cf_links[id]; // from=coarse, to=fine
            const std::array<int,6> k2{L.to_ijk[0],L.to_ijk[1],L.to_ijk[2],
                                       L.from_ijk[0],L.from_ijk[1],L.from_ijk[2]};
            const auto it = gfc.find(k2);
            if (it != gfc.end()) ++g_hit; else { ++g_miss; ++nmiss; }
            const double coeff = (it!=gfc.end() ? it->second : 0.0) / volratio;
            cf_links[id].coeff = coeff;
            new_total += coeff;
        }
        a->M.p[r0.n] += (r0.old - new_total); // diagonal: T_old -> conservative T_cf
        slot(r0.n, r0.axis, r0.high) = 0.0;

        if (cprobe && nmiss > 0)
        {
            const auto& L0 = cf_links[ids[0]];
            std::cout << "  [cfcoarse] MISS coarse (" << L0.from_ijk[0] << "," << L0.from_ijk[1]
                      << "," << L0.from_ijk[2] << ") axis=" << r0.axis << " high=" << r0.high
                      << "  fine_partners=" << ids.size() << " missing=" << nmiss
                      << "  coeff_total=" << new_total << std::endl;
        }
    }
    if (cprobe)
    {
        long tot_hit = g_hit, tot_miss = g_miss;
        MPI_Allreduce(MPI_IN_PLACE, &tot_hit,  1, MPI_LONG, MPI_SUM, pgc->mpi_comm);
        MPI_Allreduce(MPI_IN_PLACE, &tot_miss, 1, MPI_LONG, MPI_SUM, pgc->mpi_comm);
        if (p->mpirank == 0)
            std::cout << "  [cfcoarse] coarse->fine gfc lookups: hit=" << tot_hit
                      << "  MISS=" << tot_miss
                      << "  (miss => coarse row decoupled from fine -> reflux injects divergence)"
                      << std::endl;
    }
#endif
}

void hypre_ssamg::cf_velocity_correction(lexer* p, fdm* a, ghostcell* pgc,
                                         field& pcorr, field& u, field& v, field& w,
                                         double alpha)
{
#if USE_AMREX
    // REEF_CF_PROJECTION_GROUP member (2) — keep the C-F d_cf / fine-density / volume
    // convention in sync with the rest of the group; canonical note in amr_cf_coefficients.
    //
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

    // REEF_CF_CORNER_PROBE: at each coarse C-F face, compare the matrix-adjoint coarse flux
    // (Sum coeff_cf*(phiF-phiC)) against the reflux-implied flux (avg of the corrected fine face
    // velocities / (alpha*dt*DXN_c)). If the ratio is ~1 the reflux already equals the coupling
    // (a direct coarse-face correction would change nothing); a ratio != 1 is the factor the
    // reflux is off by at that face. Keyed by coarse {to_ijk, axis, high}. Per-rank partial sums.
    const bool corner = (std::getenv("REEF_CF_CORNER_PROBE") != nullptr);
    const double volratio_cp = double(rr[0])*double(rr[1])*double(rr[2]);
    std::map<std::array<int,5>, std::array<double,6>> cornacc; // {mflux, sum_ufine, count, DXN_c, sum|coeff|, sum|grad|}
    double gmin_coeff = 1e300, gmax_coeff = 0.0;   // global |L.coeff| range over fine->coarse links
    int craw = 0;                                  // raw pf/pc dump counter

    // --- TEMPORARY localization probe (REEF_CF_DIV), remove once the C-F ring is fixed ---
    // Measures the post-correction discrete divergence of each fine interface cell, split by
    // link type (high/low/both) and lateral wall distance, and reports the global worst cell.
    const bool dbg = (std::getenv("REEF_CF_DIV") != nullptr);
    const bool cf_count = (std::getenv("REEF_CF_COUNT") != nullptr);
    double g_hi = 0.0, g_lo = 0.0, g_both = 0.0;   // max|D.u| by link type (interface cells)
    double g_wall = 0.0, g_inner = 0.0;             // max|D.u| by lateral wall distance (<=1 cell)
    double l_worst = 0.0; int l_lev = -1, l_ijk[3] = {-1,-1,-1}; bool l_hi = false, l_lo = false;
    int l_fc = 0, l_fzm = 0, l_fzp = 0, l_fxm = 0, l_fxp = 0; // flag4 at worst cell + z/x neighbours
    double g_wall_noif = 0.0;                        // max|D.u| at lateral-wall NON-interface fluid cells
    double n_worst = 0.0; int n_lev = -1, n_ijk[3] = {-1,-1,-1};
    double g_cwall = 0.0;                            // coarse(lev0) lateral-wall, EXCLUDING covered cells
    double c_worst = 0.0; int c_ijk[3] = {-1,-1,-1};

    // --- Wall-cell decomposition probe (REEF_WALL_PROBE) -----------------------------
    // At the worst coarse(lev0) lateral-wall cell, capture the full velocity/pressure
    // decomposition so the divergence source is unambiguous: the 6 staggered face
    // velocities, the per-axis divergence split, this cell's pcorr and its 6 neighbour
    // pcorrs, and the 6 neighbour flag4. For a no-slip wall the wall-normal face velocity
    // must be 0; a nonzero value there is an uncorrectable source (velocity BC bug). If
    // the wall face is 0 but D.u != 0, the source is a G!=L gradient mismatch.
    const bool wallprobe = (std::getenv("REEF_WALL_PROBE") != nullptr);
    double cw_u[2] = {0,0}, cw_v[2] = {0,0}, cw_w[2] = {0,0};   // [low,high] faces
    double cw_dax[3] = {0,0,0};                                  // per-axis divergence contribution
    double cw_pc = 0.0, cw_pcn[6] = {0,0,0,0,0,0};               // cell + (x-,x+,y-,y+,z-,z+) pcorr
    int    cw_fln[6] = {0,0,0,0,0,0};                            // (x-,x+,y-,y+,z-,z+) flag4

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

                // Interface axis/side in the fine cell's own index space.
                const int axis = L.axis;
                const bool high = L.high;

                const double absnew = std::fabs(L.coeff);  // C-F face coeff over dcf

                const double pf = pf_arr(fiv);
                const double pc = pc_arr(tiv);

                if (corner && craw < 10)
                {
                    ++craw;
                    const amrex::IntVect cfiv = amrex::coarsen(fiv, rr);
                    std::cout << "  [cfraw] fiv=(" << fiv[0] << "," << fiv[1] << "," << fiv[2]
                              << ") coarsen(fiv)=(" << cfiv[0] << "," << cfiv[1] << "," << cfiv[2]
                              << ") tiv=(" << tiv[0] << "," << tiv[1] << "," << tiv[2]
                              << ") axis=" << axis << " high=" << high
                              << "  pf=" << pf << "  pc=" << pc
                              << "  pc(coarsen(fiv))=" << pc_arr(cfiv) << std::endl;
                }

                amrex::IntVect e(0, 0, 0); e[axis] = 1;

                // Sole-writer C-F correction. The corrected face is the fine cell's own face on
                // the C-F side: the face AT fiv for a high link, the face at fiv-e for a low link.
                // The interior velcorr_amrex loop skips high C-F faces (cf_masks) and never visits
                // low C-F faces (fiv-e lies outside the patch), so this applies the full
                // matrix-consistent gradient absnew*(p_coarse - p_fine) over dcf directly -- no
                // undo of an interior contribution is needed.
                const amrex::IntVect face = high ? fiv : fiv - e;
                const double grad = high ? (pc - pf) : (pf - pc);
                const double delta = alpha * p->dt * dx[axis] * absnew * grad;

                // --- TEMPORARY one-cell adjoint probe (env REEF_CF_FACE_PROBE) ---
                // Dump the raw face inputs at the corner cell (8,8,23) [all 3 of its C-F links]
                // and two non-corner patch-top cells (12,8,23),(16,8,23) [+z link], so the
                // velocity-side flux (|L.coeff|*grad) can be compared corner vs interior and
                // against the matrix Lp from projprobe. Fires both in production (solved pcorr)
                // and in projcheck (analytic pcorr -- identify by pc/pf matching projprobe).
                if (std::getenv("REEF_CF_FACE_PROBE"))
                {
                    // Target the projcheck worst cell (79,0,1) and its immediate x-neighbours,
                    // dumping EVERY C-F link that touches them (any axis/side) + the corrected face.
                    const bool tgt = (fiv[1]==0 && fiv[2]==1 && (fiv[0]>=78 && fiv[0]<=80));
                    if (tgt)
                        std::cout << "  [cffaceprobe] fiv=(" << fiv[0] << "," << fiv[1] << "," << fiv[2] << ")"
                                  << " axis=" << axis << " high=" << high
                                  << "  corrected_face=(" << face[0] << "," << face[1] << "," << face[2] << ")"
                                  << "  |L.coeff|=" << absnew << "  dx=" << dx[axis]
                                  << "  pc=" << pc << " pf=" << pf << " grad=" << grad
                                  << "  delta=" << delta
                                  << "  |L.coeff|*grad=" << (absnew * grad) << std::endl;
                }

                if      (axis == 0) u_arr(face) -= delta;
                else if (axis == 1) v_arr(face) -= delta;
                else                w_arr(face) -= delta;

                if (corner)
                {
                    const double ucorr = (axis==0)? u_arr(face) : (axis==1? v_arr(face) : w_arr(face));
                    auto& acc = cornacc[{tiv[0],tiv[1],tiv[2],axis,high?1:0}];
                    acc[0] += (absnew/volratio_cp) * (pf - pc);   // Sum coeff_cf*(phiF-phiC)
                    acc[1] += ucorr;                              // corrected fine face vel (feeds reflux)
                    acc[2] += 1.0;
                    acc[3]  = p->amrex_geometry[lev-1].CellSize(axis); // coarse cell size DXN_c
                    acc[4] += absnew;                             // |L.coeff| (matrix C-F coupling)
                    acc[5] += std::fabs(pf - pc);                 // |phiF-phiC| gradient magnitude
                    gmin_coeff = std::min(gmin_coeff, absnew);
                    gmax_coeff = std::max(gmax_coeff, absnew);
                }

                // REEF_CF_COUNT: register this C-F face write in the shared per-face counter so
                // velcorr can flag any face also written by the interior loop (count==2).
                if (cf_count) ++cf_wcount[{lev, face[0], face[1], face[2], axis}];
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
                const auto pa = pcorr.GetMultiFab(0).const_array(mfi);
                for (int ck = bx.smallEnd(2); ck <= bx.bigEnd(2); ++ck)
                for (int cj = bx.smallEnd(1); cj <= bx.bigEnd(1); ++cj)
                for (int ci = bx.smallEnd(0); ci <= bx.bigEnd(0); ++ci)
                {
                    const bool onwall = (ci==dl0[0] || ci==dh0[0]
                                       || (p->j_dir && (cj==dl0[1] || cj==dh0[1])));
                    if (!onwall) continue;
                    if (fa(ci,cj,ck) < AIR_FLAG) continue;
                    if (covered.contains(amrex::IntVect(ci,cj,ck))) continue;   // skip cells under fine patch
                    const double ddx = (ua(ci,cj,ck) - ua(ci-1,cj,ck)) / dx0[0];
                    const double ddz = (wa(ci,cj,ck) - wa(ci,cj,ck-1)) / dx0[2];
                    const double ddy = p->j_dir ? (va(ci,cj,ck) - va(ci,cj-1,ck)) / dx0[1] : 0.0;
                    const double d  = ddx + ddy + ddz;
                    const double ad = std::fabs(d);
                    if (ad > g_cwall)  g_cwall = ad;
                    if (ad > c_worst)
                    {
                        c_worst=ad; c_ijk[0]=ci; c_ijk[1]=cj; c_ijk[2]=ck;
                        if (wallprobe)
                        {
                            cw_u[0]=ua(ci-1,cj,ck); cw_u[1]=ua(ci,cj,ck);
                            cw_v[0]=va(ci,cj-1,ck); cw_v[1]=va(ci,cj,ck);
                            cw_w[0]=wa(ci,cj,ck-1); cw_w[1]=wa(ci,cj,ck);
                            cw_dax[0]=ddx; cw_dax[1]=ddy; cw_dax[2]=ddz;
                            cw_pc=pa(ci,cj,ck);
                            cw_pcn[0]=pa(ci-1,cj,ck); cw_pcn[1]=pa(ci+1,cj,ck);
                            cw_pcn[2]=pa(ci,cj-1,ck); cw_pcn[3]=pa(ci,cj+1,ck);
                            cw_pcn[4]=pa(ci,cj,ck-1); cw_pcn[5]=pa(ci,cj,ck+1);
                            cw_fln[0]=fa(ci-1,cj,ck); cw_fln[1]=fa(ci+1,cj,ck);
                            cw_fln[2]=fa(ci,cj-1,ck); cw_fln[3]=fa(ci,cj+1,ck);
                            cw_fln[4]=fa(ci,cj,ck-1); cw_fln[5]=fa(ci,cj,ck+1);
                        }
                    }
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
        {
            std::cout << "  [cf-div] worst coarse wall |D.u|=" << r_cworst
                      << " at lev=0 (" << c_ijk[0] << "," << c_ijk[1] << "," << c_ijk[2] << ")" << std::endl;
            if (wallprobe)
            {
                std::cout << "  [wallprobe] (" << c_ijk[0] << "," << c_ijk[1] << "," << c_ijk[2] << ")"
                          << "  u[lo,hi]=" << cw_u[0] << "," << cw_u[1]
                          << "  v[lo,hi]=" << cw_v[0] << "," << cw_v[1]
                          << "  w[lo,hi]=" << cw_w[0] << "," << cw_w[1] << std::endl;
                std::cout << "  [wallprobe] D.u split  dudx=" << cw_dax[0]
                          << "  dvdy=" << cw_dax[1] << "  dwdz=" << cw_dax[2] << std::endl;
                std::cout << "  [wallprobe] pcorr c=" << cw_pc
                          << "  x-/x+=" << cw_pcn[0] << "/" << cw_pcn[1]
                          << "  y-/y+=" << cw_pcn[2] << "/" << cw_pcn[3]
                          << "  z-/z+=" << cw_pcn[4] << "/" << cw_pcn[5] << std::endl;
                std::cout << "  [wallprobe] flag4  x-/x+=" << cw_fln[0] << "/" << cw_fln[1]
                          << "  y-/y+=" << cw_fln[2] << "/" << cw_fln[3]
                          << "  z-/z+=" << cw_fln[4] << "/" << cw_fln[5]
                          << "  (wall-normal face velocity must be 0 for no-slip)" << std::endl;
            }
        }
    }

    if (corner)
    {
        std::cout << "  [cfcorner] |L.coeff| range over fine->coarse links on this rank:  min="
                  << gmin_coeff << "  max=" << gmax_coeff
                  << "   (==0 => matrix has NO C-F coupling)" << std::endl;
        int shown = 0;
        for (const auto& [key, acc] : cornacc)
        {
            const double mflux    = acc[0];                        // matrix coarse coupling flux
            const double refl_vel = acc[1] / acc[2];               // reflux coarse face velocity
            const double refl_flux= refl_vel / (alpha*p->dt*acc[3]); // reflux-implied coarse flux
            const double avgcoeff = acc[4] / acc[2];               // avg |L.coeff| on this face
            const double avggrad  = acc[5] / acc[2];               // avg |phiF-phiC|
            if (shown++ < 12)
                std::cout << "  [cfcorner] coarse (" << key[0] << "," << key[1] << "," << key[2]
                          << ") axis=" << key[3] << " high=" << key[4] << " N=" << int(acc[2])
                          << "  |coeff|=" << avgcoeff << "  |grad|=" << avggrad
                          << "  matrix_flux=" << mflux << "  reflux_flux=" << refl_flux << std::endl;
        }
    }
#endif
}

void hypre_ssamg::velcorr_amr(lexer* p, fdm* a, ghostcell* pgc, density* pd,
                              field& pcorr, field& u, field& v, field& w, double alpha)
{
#if USE_AMREX
    // REEF_CF_PROJECTION_GROUP member — the complete multi-level projection velocity
    // correction, owned by the solver so it reads the assembled operator directly.
    //
    // Interior faces are sourced from a->M (the exact adjoint of the divergence the solve
    // saw): |M.slot|*DXN == CPOR*POR/(roface*DXP), bit-identical to the old ucorr/vcorr/wcorr
    // on fluid faces, but auto-zero exactly where amr_cf_coefficients zeroed the slot (C-F
    // seams, high AND low, at every corner) and where poisson_pcorr folded a wall face. That
    // removes the fragile cf_masks: the skip-set IS the zero-set. a->Mrow (persistent) maps
    // each solved cell to its row n; a->M is read live (nothing between the solve and this
    // call rewrites it — reference_start only shifts the pressure gauge).
    //
    // Free surface: poisson moved the air-face coupling to the RHS (M.slot==0) but the
    // Dirichlet correction against the pcorr ghost is real, so a neighbour == AIR_FLAG face
    // falls back to the roface gradient form (identical to the single-level ucorr). This is
    // the one branch that needs the model face density, hence pd. Domain walls that carry the
    // AIR_FLAG ghost also take this branch and are then overwritten by start1/2/3, exactly as
    // the old ucorr behaved.
    LOOP
    {
        const int n = a->Mrow(i,j,k);
        if(n < 0) continue;

        // x+ (high) face -> u(i,j,k), couples cells (i,j,k) and (i+1,j,k)
        {
            const double dp = pcorr(i+1,j,k) - pcorr(i,j,k);
            if(a->M.n[n] != 0.0)
                u(i,j,k) += alpha*p->dt * a->M.n[n] * p->DXN[IP] * dp;
            else if(p->flag4(i+1,j,k)==AIR_FLAG)
                u(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(dp/(p->DXP[IP]*a->rofx(i,j,k)));
        }

        // y+ (high) face -> v(i,j,k)
        if(p->j_dir==1)
        {
            const double dp = pcorr(i,j+1,k) - pcorr(i,j,k);
            if(a->M.w[n] != 0.0)
                v(i,j,k) += alpha*p->dt * a->M.w[n] * p->DYN[JP] * dp;
            else if(p->flag4(i,j+1,k)==AIR_FLAG)
                v(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(dp/(p->DYP[JP]*a->rofy(i,j,k)));
        }

        // z+ (high) face -> w(i,j,k)
        {
            const double dp = pcorr(i,j,k+1) - pcorr(i,j,k);
            if(a->M.t[n] != 0.0)
                w(i,j,k) += alpha*p->dt * a->M.t[n] * p->DZN[KP] * dp;
            else if(p->flag4(i,j,k+1)==AIR_FLAG)
                w(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*(dp/(p->DZP[KP]*a->rofz(i,j,k)));
        }
    }

    // Coarse-fine interface (fine sub-faces, gradient/G side) then reflux the covered
    // coarse faces to the area-summed fine flux: together D(1/rho)G = L at the interface.
    cf_velocity_correction(p,a,pgc,pcorr,u,v,w,alpha);
    cf_average_down_velocity(p,u,v,w);
#endif
}

void hypre_ssamg::cf_average_down_velocity(lexer* p, field& u, field& v, field& w)
{
#if USE_AMREX
    // Sync the coarse velocity under each fine patch with the fine solution: build face
    // MultiFabs from the staggered velocity, average_down_faces (fine -> covered coarse faces,
    // incl. the C-F interface face), copy the synced coarse faces back. This makes the coarse
    // C-F face velocity equal the sum/average of the fine sub-fluxes (reflux) and the covered
    // coarse cells equal the fine average. (Kept in sync with the pjm_corr predictor-side copy.)
    const int jdir = p->j_dir;
    field* vfld[AMREX_SPACEDIM] = {&u, &v, &w};

    for (int lev = p->nlevs-1; lev >= 1; --lev)
    {
        const int clev = lev-1;
        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> ffine, fcrse;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            ffine[dir].define(amrex::convert(p->amrex_box_array[lev],  e), p->amrex_distribution_mapping[lev],  1, 0);
            fcrse[dir].define(amrex::convert(p->amrex_box_array[clev], e), p->amrex_distribution_mapping[clev], 1, 0);

            // staggered cell-centred velocity -> face MultiFab (face f = vel(f-e))
            auto fill_face = [&] (amrex::MultiFab& fmf, int amrlev)
            {
                auto& v_mf = vfld[dir]->GetMultiFab(amrlev);
                for (amrex::MFIter mfi(fmf); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& fbx = mfi.validbox();
                    auto       f = fmf.array(mfi);
                    const auto vv = v_mf.const_array(mfi);
                    amrex::LoopOnCpu(fbx, [&] (int i, int j, int k) noexcept
                    { f(i,j,k) = vv(i-e[0], j-e[1], k-e[2]); });
                }
            };
            fill_face(ffine[dir], lev);
            fill_face(fcrse[dir], clev);   // pre-fill so uncovered coarse faces survive
        }

        const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> ffp{&ffine[0], &ffine[1], &ffine[2]};
        const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>       fcp{&fcrse[0], &fcrse[1], &fcrse[2]};

        amrex::average_down_faces(ffp, fcp, p->ref_vec, p->amrex_geometry[clev]);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if (dir==1 && jdir!=1) continue;
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            auto& v_mf = vfld[dir]->GetMultiFab(clev);
            for (amrex::MFIter mfi(v_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                auto       vv = v_mf.array(mfi);
                const auto f  = fcrse[dir].const_array(mfi);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
                { vv(i,j,k) = f(i+e[0], j+e[1], k+e[2]); });
            }
        }
    }
#endif
}

void hypre_ssamg::cf_lowface_save_restore(lexer* p, field& u, field& v, field& w, bool save)
{
#if USE_AMREX
    // The fine cell's LOW C-F normal face is the ghost face (fiv-e). cf_velocity_correction writes
    // the matrix-consistent gradient there; a subsequent start1/2/3 (FillPatchTwoLevels +
    // FillCoarseFineNormalGhost) overwrites it with coarse interpolation, breaking D.u = -L.pcorr
    // at the fine cell. Save these faces before such a start and restore them after, so the
    // correction survives into the divergence/advection that consumes it.
    if (p->nlevs <= 1) return;

    if (save) cf_lowface_store.clear();

    field* vfld[3] = {&u, &v, &w};

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        for (const auto& L : cf_links)
        {
            if (L.from_part != lev)         continue; // fine cell on this level
            if (L.to_part   >= L.from_part) continue; // fine->coarse couplings only
            if (L.high)                     continue; // LOW faces only (HIGH are valid fine faces)

            const int axis = L.axis;
            amrex::IntVect e(0,0,0); e[axis] = 1;
            const amrex::IntVect face(L.from_ijk[0]-e[0], L.from_ijk[1]-e[1], L.from_ijk[2]-e[2]);

            auto& mf = vfld[axis]->GetMultiFab(lev);
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                if (!mfi.fabbox().contains(face)) continue;
                auto arr = mf.array(mfi);
                const std::array<int,5> key{lev, face[0], face[1], face[2], axis};
                if (save)
                    cf_lowface_store[key] = arr(face);
                else
                {
                    const auto it = cf_lowface_store.find(key);
                    if (it != cf_lowface_store.end()) arr(face) = it->second;
                }
                break; // face is owned by one box on this level
            }
        }
    }
#endif
}

void hypre_ssamg::cf_velocity_fill_from_coarse(lexer* p, fdm* a, ghostcell* pgc,
                                               field& u, field& v, field& w)
{
#if USE_AMREX
    // Slave each fine-level HIGH coarse-fine normal face to the coarse value beneath it:
    //   u_fine(fiv) = u_coarse(coarsen(fiv)).
    // The HIGH C-F face is a *valid* fine face the predictor freely updates -- exactly where the
    // hydrostatic band residual leaks. The coarse C-F face is interior to the coarse level (no
    // leak), so the fine face inherits the clean, balanced coarse velocity instead of the
    // predictor's leaked value, and the rhs divergence the solver sees is clean. (LOW C-F faces
    // are fine ghosts already filled from coarse by start1, so they need no action here.)
    // Injection keeps avg(fine sub-faces) == coarse, so it is consistent with the fine->coarse
    // reflux in cf_average_down_velocity -- and MUST run BEFORE that reflux, so the coarse source
    // is the coarse predictor value, not the averaged-down (leaked) one.
    if (p->nlevs <= 1) return;

    const auto& rr = p->ref_vec;
    field* vfld[3] = {&u, &v, &w};

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const amrex::BoxArray cfba = amrex::coarsen(u.GetMultiFab(lev).boxArray(), rr);

        // Coarse normal velocity on each fine box's coarsened layout (read as-is, BEFORE any
        // average-down, so it is the clean coarse predictor value).
        amrex::MultiFab cvel[3];
        for (int d = 0; d < 3; ++d)
        {
            if (d == 1 && p->j_dir != 1) continue;
            cvel[d].define(cfba, u.GetMultiFab(lev).DistributionMap(), 1, 1);
            cvel[d].setVal(0.0);
            cvel[d].ParallelCopy(vfld[d]->GetMultiFab(lev - 1), 0, 0, 1,
                                 amrex::IntVect(0), amrex::IntVect(1),
                                 p->amrex_geometry[lev - 1].periodicity());
        }

        auto& u_mf = u.GetMultiFab(lev);
        auto& v_mf = v.GetMultiFab(lev);
        auto& w_mf = w.GetMultiFab(lev);

        for (amrex::MFIter mfi(u_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            auto u_arr = u_mf.array(mfi);
            auto v_arr = v_mf.array(mfi);
            auto w_arr = w_mf.array(mfi);
            const auto cvu = cvel[0].const_array(mfi);
            const auto cvw = cvel[2].const_array(mfi);
            amrex::Array4<const double> cvv;
            if (p->j_dir == 1) cvv = cvel[1].const_array(mfi);

            for (const auto& L : cf_links)
            {
                if (L.from_part != lev || L.to_part >= L.from_part) continue; // fine->coarse only
                if (!L.high) continue;                                        // valid fine faces only
                const int axis = L.axis;
                if (axis == 1 && p->j_dir != 1) continue;

                const amrex::IntVect fiv(L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]);
                if (!bx.contains(fiv)) continue;

                // covered coarse cell containing fiv; u_coarse there IS the coarse C-F face.
                const amrex::IntVect cf = amrex::coarsen(fiv, rr);

                if      (axis == 0) u_arr(fiv) = cvu(cf);
                else if (axis == 2) w_arr(fiv) = cvw(cf);
                else                v_arr(fiv) = cvv(cf);
            }
        }
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
