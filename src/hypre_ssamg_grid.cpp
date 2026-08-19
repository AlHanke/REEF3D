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
#include <array>
#include <map>
#include <set>

void hypre_ssamg::make_grid_7p(lexer *p, fdm *a, ghostcell *pgc)
{
    // Grid setup
    #if USE_AMREX
    numparts   = p->nlevs;
    #else
    numparts   = 1;
    #endif
    dimensions = 3;
    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);

    // Variables
    numvar = 1;
    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;

    for (int lev = 0; lev < numparts; ++lev)
        HYPRE_SStructGridSetVariables(grid, lev, numvar, vartypes);


    // Extends
    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;

    iupper[0] = p->knox + p->origin_i - 1;
    iupper[1] = p->knoy + p->origin_j - 1;
    iupper[2] = p->knoz + p->origin_k - 1;

    HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);

    #if USE_AMREX
    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        for (int b = 0; b < (int)p->amrex_box_array[lev].size(); ++b)
        {
            if (p->amrex_distribution_mapping[lev][b] != p->mpirank) continue;
            const auto& box = p->amrex_box_array[lev][b];
            int lo[3] = {box.smallEnd(0), box.smallEnd(1), box.smallEnd(2)};
            int hi[3] = {box.bigEnd(0),   box.bigEnd(1),   box.bigEnd(2)};
            HYPRE_SStructGridSetExtents(grid, lev, lo, hi);
        }
    }
    #endif

    // Periodicity -- must be set on EVERY part before Assemble.
    //
    // REEF3D handles periodic BCs by exchanging ghost cells, so poisson_pcorr leaves the
    // off-diagonal coefficient intact on a cell at a periodic boundary (its neighbour is a
    // valid fluid ghost, not a wall, so none of the fold/Dirichlet branches fire). HYPRE
    // therefore receives a stencil entry pointing outside the part's index space and, without
    // this call, silently DROPS it -- the wrap-around coupling vanishes from the operator and
    // the periodic direction is solved as if it had a free (zero-Neumann-ish) boundary.
    // hypre_struct/hypre_struct_amr already do this; the SStruct path did not.
    //
    // The period is the part's OWN index-space length: level 0 spans the base domain, an AMR
    // level lev spans it refined by ref^lev, so each part gets its own value.
    //
    // Setting it unconditionally on every part is safe, including for interior fine patches.
    // Periodicity here is a statement about the part's index space, not a claim that the patch
    // sits on a physical periodic boundary: hypre resolves an out-of-part stencil entry by
    // shifting it by +/- the period and looking for a box there. Two cases:
    //   * the patch does NOT reach both ends of the periodic direction -- the shifted index lands
    //     entirely outside the part's boxes, hypre finds nothing and drops the entry exactly as it
    //     does today. Inert; no spurious wrap-around can be manufactured, because the shift is the
    //     FULL refined domain length and only a box near the opposite extreme can match it.
    //   * the patch DOES span the periodic direction (the common case -- a refined band following
    //     the free surface across a periodic channel) -- the low- and high-edge boxes are wired to
    //     each other, which is the correct coupling and the whole point of the call.
    // The one gap is a patch that touches ONE periodic end only: its wrap-around neighbour is a
    // COARSE cell, i.e. a C-F coupling, and the cf_links halo construction filters candidates with
    // in_domain(), which rejects a halo sitting outside the level's domain box. So that coupling is
    // still missed. That is a pre-existing limitation of the C-F link construction (periodic C-F
    // interfaces are not built at all), not something this call introduces or worsens.
    //
    // NOTE (HYPRE docs): some hypre solvers impose power-of-two restrictions on the size of a
    // periodic dimension. SSAMG's semi-coarsening is one -- a period that stops being even limits
    // how far the periodic direction coarsens. This only bites the single-level SSAMG path; the
    // nlevs>1 path assembles to ParCSR and coarsens algebraically in BoomerAMG, where the grid
    // periodicity only shapes the assembled matrix graph.
    {
        HYPRE_Int periodic[3];

        for (int lev = 0; lev < numparts; ++lev)
        {
            #if USE_AMREX
            const amrex::IntArray amrex_periodic = p->amrex_geometry[lev].isPeriodic();
            const amrex::Box &dom = p->amrex_geometry[lev].Domain();
            periodic[0] = amrex_periodic[0] ? dom.length(0) : 0;
            periodic[1] = amrex_periodic[1] ? dom.length(1) : 0;
            periodic[2] = amrex_periodic[2] ? dom.length(2) : 0;
            #else
            periodic[0] = p->periodic1 > 0 ? p->gknox : 0;
            periodic[1] = p->periodic2 > 0 ? p->gknoy : 0;
            periodic[2] = p->periodic3 > 0 ? p->gknoz : 0;
            #endif

            HYPRE_SStructGridSetPeriodic(grid, lev, periodic);
        }
    }

    HYPRE_SStructGridAssemble(grid);


    // 7-point stencil: centre + ±x ±y ±z
    int offsets[stencil_size][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};
    variable = 0;

    HYPRE_SStructStencilCreate(dimensions, stencil_size, &stencil);
    for (int entry = 0; entry < 7; ++entry)
        HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);


    // Single level -> SSAMG (native SStruct). Multi level -> assemble as ParCSR so the
    // graph-coupled multi-part operator can be solved with PCG+BoomerAMG (SSAMG
    // setup truncates on multi-part grids).
    #if USE_AMREX
    object_type = (p->nlevs > 1) ? HYPRE_PARCSR : HYPRE_SSTRUCT;
    #else
    object_type = HYPRE_SSTRUCT;
    #endif

    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetObjectType(graph, object_type);
    for (int lev = 0; lev < numparts; ++lev)
        HYPRE_SStructGraphSetStencil(graph, lev, variable, stencil);

    #if USE_AMREX
    // State shared across all levels: dedup guard, per-cell non-stencil counter and
    // the coupling records. A cell on an intermediate part can receive entries in two
    // level iterations (as fine cell toward the coarser level, then as coarse cell
    // toward the finer one), so these must persist so the entry indices stay contiguous.
    cf_links.clear();
    cf_masks.clear();
    std::set<std::array<int,stencil_size>>     added;   // (from_part, from_ijk, to_ijk)
    std::map<std::array<int,4>,int> nnz;     // (from_part, from_ijk) -> #non-stencil so far

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const auto& fine_ba = p->amrex_box_array[lev];
        const auto& rr      = p->ref_vec;

        auto base_ba = p->amrex_box_array[0];
        for (int n = 0; n < lev; ++n) base_ba.refine(rr);
        const auto base_hi = amrex::IntVect(amrex::ubound(base_ba.minimalBox()));

        auto in_domain = [&](const amrex::Box& b)
        {
            for (int i = 0; i < 3; ++i)
                if (b.bigEnd(i) < 0 || b.smallEnd(i) > base_hi[i]) return false;

            return true;
        };

        std::vector<std::tuple<amrex::Box, amrex::Box, int>> cf_pairs;

        for (int b = 0; b < (int)fine_ba.size(); ++b)
        {
            for (auto halo : amrex::boxDiff(amrex::grow(fine_ba[b], 1), fine_ba[b]))
            {
                if (!in_domain(halo)) continue;

                // Identify the face normal as the single thin direction. In pseudo-2D the
                // inactive y-direction is a singleton cell, so an x- or z-normal halo is thin
                // in BOTH its normal and y -> n_thin would be 2 and the face (e.g. the x-pos
                // C-F interface) would be wrongly skipped. Exclude the inactive y-axis from the
                // count so pseudo-2D x/z faces are detected as thin-in-one.
                int normal = -1, n_thin = 0;
                for (int i = 0; i < 3; ++i)
                {
                    if (i == 1 && !p->j_dir) continue;   // y inactive (pseudo-2D): not a face normal
                    if (halo.smallEnd(i) == halo.bigEnd(i)) { normal = i; ++n_thin; }
                }
                if (n_thin != 1) continue;

                halo &= amrex::grow(fine_ba[b], normal, 1);

                // A genuine coarse-fine ghost cell lies over coarse (uncovered) territory. Where
                // this fine box's halo shell overlaps ANOTHER fine box at the same level the seam
                // is fine-fine, not coarse-fine: HYPRE's intra-part 7-point stencil already
                // couples the two fine boxes there, so emitting a C-F coupling would zero that
                // stencil slot (amr_cf_coefficients) and decouple the boxes across the seam,
                // injecting a spurious pressure jump / stray interface velocity. Restrict the halo
                // to the part not covered by fine_ba so only real C-F faces produce cf_pairs.
                for (const amrex::Box& halo_cf : fine_ba.complementIn(halo))
                {
                    const auto box_c = amrex::coarsen(halo_cf, rr);
                    const auto box_f = amrex::refine(box_c, rr);

                    // Emit ONLY along the halo's normal axis. box_c = coarsen(halo_cf) is a single
                    // coarse layer in `normal` (the halo is thin there -- n_thin==1), which is the
                    // assumption `cn = coarse_box.smallEnd(normal)` (line ~309) relies on. Iterating
                    // the transverse axes i!=normal paired box_c -- which spans MANY cells in i --
                    // with fine cells via smallEnd(i)=0, coupling interface fine cells to a coarse
                    // cell at the domain edge (the z=0 / x=0 corruption that broke both the matrix
                    // C-F coupling and the velocity corrector -> spurious C-F velocity -> geyser).
                    // A staircase interface's other-axis C-F faces are caught by their own
                    // thin-in-that-axis halos, so no real coupling is lost.
                    const int i = normal;
                    for (auto& [idx, overlap] : fine_ba.intersections(amrex::adjCellHi(box_f, i, 1)))
                        cf_pairs.emplace_back(box_c, overlap, 2*i+1);
                    for (auto& [idx, overlap] : fine_ba.intersections(amrex::adjCellLo(box_f, i, 1)))
                        cf_pairs.emplace_back(box_c, overlap, 2*i);
                }
            }
        }

        // REEF_GRID_DEBUG: dump the fine box layout, the halo faces each box produces, and the
        // cf_pairs dir-histogram, so a missing C-F face (e.g. the x-pos edge yielding no links)
        // can be traced to the halo/in_domain/complementIn step. Rank 0 only.
        if (std::getenv("REEF_GRID_DEBUG") && p->mpirank == 0)
        {
            std::cout << "\n  [griddbg] lev=" << lev << " base_hi=("
                      << base_hi[0] << "," << base_hi[1] << "," << base_hi[2] << ")"
                      << "  nfineboxes=" << fine_ba.size() << std::endl;
            for (int b = 0; b < (int)fine_ba.size(); ++b)
            {
                const auto& fb = fine_ba[b];
                std::cout << "  [griddbg] box " << b << " ["
                          << fb.smallEnd(0) << ".." << fb.bigEnd(0) << " , "
                          << fb.smallEnd(1) << ".." << fb.bigEnd(1) << " , "
                          << fb.smallEnd(2) << ".." << fb.bigEnd(2) << "]" << std::endl;
                for (auto halo : amrex::boxDiff(amrex::grow(fb, 1), fb))
                {
                    const bool ind = in_domain(halo);
                    int normal=-1, n_thin=0;
                    for (int i=0;i<3;++i)
                    {
                        if (i==1 && !p->j_dir) continue;
                        if (halo.smallEnd(i)==halo.bigEnd(i)) { normal=i; ++n_thin; }
                    }
                    int ncomp = 0;
                    if (ind && n_thin==1)
                    {
                        amrex::Box h2 = halo & amrex::grow(fb, normal, 1);
                        ncomp = (int)fine_ba.complementIn(h2).size();
                    }
                    std::cout << "  [griddbg]   halo ["
                              << halo.smallEnd(0) << ".." << halo.bigEnd(0) << " , "
                              << halo.smallEnd(1) << ".." << halo.bigEnd(1) << " , "
                              << halo.smallEnd(2) << ".." << halo.bigEnd(2) << "]"
                              << " in_domain=" << ind << " normal=" << normal
                              << " n_thin=" << n_thin << " complementIn_pieces=" << ncomp << std::endl;
                }
            }
            int dirh[6]={0,0,0,0,0,0};
            for (const auto& cp : cf_pairs) { int d=std::get<2>(cp); if(d>=0&&d<6) ++dirh[d]; }
            std::cout << "  [griddbg] cf_pairs by dir (x-,x+,y-,y+,z-,z+)= "
                      << dirh[0] << " " << dirh[1] << " " << dirh[2] << " "
                      << dirh[3] << " " << dirh[4] << " " << dirh[5] << std::endl;
        }

        // --- Add the coarse-fine couplings as non-stencil graph entries ---------
        // HYPRE requires the entry to be added by the process that owns the "from"
        // cell (part,index). Levels live on independent distribution mappings, so
        // the fine cell and its coarse neighbour may be owned by different ranks;
        // ownership is therefore checked per direction. A dedup set guards against
        // the same coupling being emitted by overlapping cf_pairs (which would
        // otherwise double-count the connection in the matrix).

        auto owner_rank = [&](int level, const int iv[3]) -> int
        {
            const amrex::IntVect v(iv[0], iv[1], iv[2]);
            const auto isects = p->amrex_box_array[level].intersections(amrex::Box(v, v));
            if (isects.empty()) return -1;
            return p->amrex_distribution_mapping[level][isects[0].first];
        };

        auto add_entry = [&](int from_part, const int from[3],
                             int to_part,   const int to[3])
        {
            if (owner_rank(from_part, from) != p->mpirank) return;

            const std::array<int,7> key{from_part, from[0], from[1], from[2],
                                                   to[0],   to[1],   to[2]};
            if (!added.insert(key).second) return;

            int f[3] = {from[0], from[1], from[2]};
            int t[3] = {to[0],   to[1],   to[2]};
            HYPRE_SStructGraphAddEntries(graph, from_part, f, variable,
                                                to_part,   t, variable);

            int  axis = 0;
            bool high = false;
            amrex::IntVect fiv(from[0], from[1], from[2]);
            amrex::IntVect tiv(to[0], to[1], to[2]);
            if (from_part > to_part)
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

            // Non-stencil entry index for this from-cell: stencil_size + k-th entry.
            const std::array<int,4> ck{from_part, from[0], from[1], from[2]};
            const int k = nnz[ck]++;
            cf_links.push_back({from_part, {f[0], f[1], f[2]},
                                to_part,   {t[0], t[1], t[2]}, axis, high, stencil_size + k, 0.0});

            // Mask the fine C-F face out of the interior velcorr loop so cf_velocity_correction
            // is its SOLE writer. HIGH fine->coarse link: the fine cell's own high face (cell
            // `from`, +axis). LOW fine->coarse link: the fine cell's low face = cell (from-e)'s
            // +axis face. The low case matters on a phi-following STAIRCASE patch, where from-e is
            // a valid in-patch fine cell the interior loop would otherwise ALSO write (count==2
            // double-correction -> D.G != L). Harmless at a flat edge (from-e is out-of-patch, so
            // the interior loop never visits it). Only fine->coarse links (from_part>to_part) are
            // corrected by cf_velocity_correction, so only those are masked.
            if(from_part > to_part)
            {
                if(high)
                    cf_masks.insert({from_part, from[0], from[1], from[2], axis});
                else
                {
                    int e[3] = {0,0,0}; e[axis] = 1;
                    cf_masks.insert({from_part, from[0]-e[0], from[1]-e[1], from[2]-e[2], axis});
                }
            }
        };

        for (const auto& [coarse_box, fine_box, dir] : cf_pairs)
        {
            // dir encodes the coarse->fine direction as 2*axis + side, matching the
            // 2*i+1 (high) / 2*i (low) tags set when cf_pairs is built.
            const int  normal  = dir / 2;
            const bool fine_hi = (dir & 1);

            const int cn = coarse_box.smallEnd(normal); // single coarse layer in normal dir

            // fine -> coarse: every fine cell couples to its one coarse neighbour
            for (int k = fine_box.smallEnd(2); k <= fine_box.bigEnd(2); ++k)
            for (int j = fine_box.smallEnd(1); j <= fine_box.bigEnd(1); ++j)
            for (int i = fine_box.smallEnd(0); i <= fine_box.bigEnd(0); ++i)
            {
                const int f[3] = {i, j, k};
                const amrex::IntVect cc = amrex::coarsen(amrex::IntVect(i, j, k), rr);
                int c[3] = {cc[0], cc[1], cc[2]};
                c[normal] = cn; // neighbour lives on the coarse side of the interface

                // VALIDATION (REEF_CF_LINK_CHECK): the coarse neighbour must differ from
                // coarsen(fine) by exactly +-1 in the normal only. A larger jump (e.g. cn=0 while
                // coarsen(fine)[normal]~51) means this cf_pair matched the wrong coarse layer ->
                // the fine cell is coupled to a bogus coarse cell (broken matrix + corrector).
                if (std::getenv("REEF_CF_LINK_CHECK"))
                {
                    const int dn = c[normal] - cc[normal];
                    if (dn != 1 && dn != -1)
                        std::cout << "  [cflinkchk] BAD pair: fine (" << i << "," << j << "," << k
                                  << ") coarsen=(" << cc[0] << "," << cc[1] << "," << cc[2]
                                  << ") -> coarse (" << c[0] << "," << c[1] << "," << c[2]
                                  << ")  normal=" << normal << " cn=" << cn << " dir=" << dir
                                  << "  coarse_box[" << coarse_box.smallEnd(normal) << ".."
                                  << coarse_box.bigEnd(normal) << "] fine_box_n["
                                  << fine_box.smallEnd(normal) << ".." << fine_box.bigEnd(normal)
                                  << "]" << std::endl;
                }

                add_entry(lev, f, lev - 1, c);
            }

            // coarse -> fine: every coarse cell couples to all fine cells across the face
            for (int k = coarse_box.smallEnd(2); k <= coarse_box.bigEnd(2); ++k)
            for (int j = coarse_box.smallEnd(1); j <= coarse_box.bigEnd(1); ++j)
            for (int i = coarse_box.smallEnd(0); i <= coarse_box.bigEnd(0); ++i)
            {
                const int c[3] = {i, j, k};
                const amrex::Box foot(amrex::IntVect(i, j, k), amrex::IntVect(i, j, k));
                amrex::Box face = fine_hi ? amrex::adjCellHi(amrex::refine(foot, rr), normal, 1)
                                          : amrex::adjCellLo(amrex::refine(foot, rr), normal, 1);
                face &= fine_box;
                if (face.isEmpty()) continue;

                for (int fk = face.smallEnd(2); fk <= face.bigEnd(2); ++fk)
                for (int fj = face.smallEnd(1); fj <= face.bigEnd(1); ++fj)
                for (int fi = face.smallEnd(0); fi <= face.bigEnd(0); ++fi)
                {
                    const int f[3] = {fi, fj, fk};
                    add_entry(lev - 1, c, lev, f);
                }
            }
        }
    }
    #endif

    HYPRE_SStructGraphAssemble(graph);

    HYPRE_SStructMatrixCreate(pgc->mpi_comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, object_type);
    HYPRE_SStructMatrixInitialize(A);

    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &b);
    HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &x);
    HYPRE_SStructVectorSetObjectType(b, object_type);
    HYPRE_SStructVectorSetObjectType(x, object_type);
    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);
}

void hypre_ssamg::destroy_grid()
{
    HYPRE_SStructStencilDestroy(stencil);
    HYPRE_SStructGraphDestroy(graph);
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructVectorDestroy(x);
}