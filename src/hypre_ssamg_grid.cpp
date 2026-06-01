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

void hypre_ssamg::make_grid_7p(lexer *p, fdm *a, ghostcell *pgc)
{
    numparts   = p->nlevs;
    part       = 0;
    dimensions = p->j_dir ? 3 : 2;
    variable   = 0;
    numvar     = 1;
    object_type = HYPRE_SSTRUCT;

    ilower[0] = p->origin_i;
    ilower[1] = p->origin_j;
    ilower[2] = p->origin_k;

    iupper[0] = p->knox + p->origin_i - 1;
    iupper[1] = p->knoy + p->origin_j - 1;
    iupper[2] = p->knoz + p->origin_k - 1;

    vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;

    HYPRE_SStructGridCreate(pgc->mpi_comm, dimensions, numparts, &grid);

    HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);
    HYPRE_SStructGridSetVariables(grid, 0, numvar, vartypes);

#if USE_AMREX
    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        HYPRE_SStructGridSetVariables(grid, lev, numvar, vartypes);
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

    HYPRE_SStructGridAssemble(grid);

    // 7-point stencil: centre + ±x ±y ±z
    int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};

    HYPRE_SStructStencilCreate(dimensions, 7, &stencil);
    for (int entry = 0; entry < 7; ++entry)
        HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], variable);

    HYPRE_SStructGraphCreate(pgc->mpi_comm, grid, &graph);
    HYPRE_SStructGraphSetObjectType(graph, object_type);
    for (int lev = 0; lev < numparts; ++lev)
        HYPRE_SStructGraphSetStencil(graph, lev, variable, stencil);
#if USE_AMREX
    amr_graph_entries(p, pgc);
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

// Explicit inter-part graph entries for the coarse-fine interface (the literal
// Figure 16 approach).  For each coarse cell on the boundary of the covered
// region we add a graph entry coupling it to each fine cell it shares a face
// with, and vice versa.  This is only necessary when the AMRNEW SetAMRPart
// machinery does not generate the inter-part connectivity automatically.
void hypre_ssamg::amr_graph_entries(lexer* p, ghostcell* pgc)
{
#if USE_AMREX
    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const int coarse_part = lev - 1;
        const int fine_part   = lev;
        const int rf0 = p->ref_vec[0];
        const int rf1 = p->ref_vec[1];
        const int rf2 = p->ref_vec[2];
        const amrex::IntVect rf_iv(rf0, rf1, rf2);

        struct FaceDir { int dir; int side; };
        const FaceDir fdirs[] = {{0,0},{0,1},{1,0},{1,1},{2,0},{2,1}};
        const int nfd = 2 * dimensions;

        // Returns true if coarse cell (ci,cj,ck) is owned by this rank.
        // Level 0 has one box per rank stored in ilower/iupper.
        // Level > 0 uses amrex_box_array / amrex_distribution_mapping.
        auto coarse_is_local = [&](int ci, int cj, int ck) -> bool
        {
            if (coarse_part == 0)
                return (ci >= ilower[0] && ci <= iupper[0] &&
                        cj >= ilower[1] && cj <= iupper[1] &&
                        ck >= ilower[2] && ck <= iupper[2]);
            for (int bb = 0; bb < (int)p->amrex_box_array[coarse_part].size(); ++bb)
            {
                if (p->amrex_distribution_mapping[coarse_part][bb] != p->mpirank) continue;
                if (p->amrex_box_array[coarse_part][bb].contains(amrex::IntVect(ci,cj,ck)))
                    return true;
            }
            return false;
        };

        for (int b = 0; b < (int)p->amrex_box_array[lev].size(); ++b)
        {
            if (p->amrex_distribution_mapping[lev][b] != p->mpirank) continue;
            const auto& fine_box = p->amrex_box_array[lev][b];
            const amrex::Box covered = amrex::coarsen(fine_box, rf_iv);

            for (int fd = 0; fd < nfd; ++fd)
            {
                const int d    = fdirs[fd].dir;
                const int side = fdirs[fd].side;
                const int d1   = (d == 0) ? 1 : 0;
                const int d2   = (d == 2) ? 1 : 2;
                // Refinement factors in the two tangential directions.
                const int rf1d = (d == 0) ? rf1 : rf0;
                const int rf2d = (d == 2) ? rf1 : rf2;

                const int coarse_d    = (side == 1) ? covered.bigEnd(d) + 1 : covered.smallEnd(d) - 1;
                const int fine_face_d = (side == 1) ? fine_box.bigEnd(d)    : fine_box.smallEnd(d);

                // Skip faces that adjoin a physical domain boundary: coarse_d
                // would fall outside the global coarse grid and no box registers it.
                const amrex::Box& coarse_domain = p->amrex_geometry[coarse_part].Domain();
                if (coarse_d < coarse_domain.smallEnd(d) || coarse_d > coarse_domain.bigEnd(d))
                    continue;

                // --- fine -> coarse ---
                // Added by the rank that owns the fine cell (this rank).
                // coarse tangential index = fine index / refinement factor.
                for (int a = fine_box.smallEnd(d1); a <= fine_box.bigEnd(d1); ++a)
                for (int c = fine_box.smallEnd(d2); c <= fine_box.bigEnd(d2); ++c)
                {
                    int fine_idx[3], coarse_idx[3];
                    fine_idx[d]    = fine_face_d;
                    fine_idx[d1]   = a;
                    fine_idx[d2]   = c;
                    coarse_idx[d]  = coarse_d;
                    coarse_idx[d1] = a / rf1d;
                    coarse_idx[d2] = c / rf2d;

                    HYPRE_SStructGraphAddEntries(graph,
                        fine_part,   fine_idx,   variable,
                        coarse_part, coarse_idx, variable);
                }

                // --- coarse -> fine ---
                // Added only by the rank that owns the coarse cell.
                // For each coarse cell on the boundary face of the covered region,
                // add entries to all fine cells in its tangential footprint.
                for (int a = covered.smallEnd(d1); a <= covered.bigEnd(d1); ++a)
                for (int c = covered.smallEnd(d2); c <= covered.bigEnd(d2); ++c)
                {
                    int coarse_idx[3];
                    coarse_idx[d]  = coarse_d;
                    coarse_idx[d1] = a;
                    coarse_idx[d2] = c;

                    if (!coarse_is_local(coarse_idx[0], coarse_idx[1], coarse_idx[2]))
                        continue;

                    const int fi_lo = a * rf1d,      fi_hi = (a + 1) * rf1d - 1;
                    const int fj_lo = c * rf2d,      fj_hi = (c + 1) * rf2d - 1;

                    for (int fa = fi_lo; fa <= fi_hi; ++fa)
                    {
                        if (fa < fine_box.smallEnd(d1) || fa > fine_box.bigEnd(d1)) continue;
                        for (int fc = fj_lo; fc <= fj_hi; ++fc)
                        {
                            if (fc < fine_box.smallEnd(d2) || fc > fine_box.bigEnd(d2)) continue;

                            int fine_idx[3];
                            fine_idx[d]  = fine_face_d;
                            fine_idx[d1] = fa;
                            fine_idx[d2] = fc;

                            HYPRE_SStructGraphAddEntries(graph,
                                coarse_part, coarse_idx, variable,
                                fine_part,   fine_idx,   variable);
                        }
                    }
                }
            }
        }
    }
#endif
}
