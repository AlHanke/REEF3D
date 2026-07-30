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

#include"ghostcell.h"
#include"lexer.h"

// ---------------------------------------------------------------------------
// gcb4_generate — build the cell-centred ghost-cell boundary list from the flag
// field and the bcside settings, instead of reading it from the DIVEMesh grid
// file.
//
// Why generate rather than read: the grid-file list (read_grid.cpp, gcb4[i][0..2]
// = isurf/jsurf/ksurf) is expressed in the LEGACY rank-local index space, where
// the global index is i + origin_i. Under AMReX the decomposition is a BoxArray
// with its own ownership, and a cell in this rank's legacy subdomain may belong
// to a different rank's box entirely — so the file list can be neither addressed
// nor tile-tagged without an MPI redistribution. It also never survives a
// regrid, which re-chops the domain underneath it.
//
// Generating from flag4 sidesteps all of that: every entry is produced inside a
// TILE_LOOP, so it is owned by the rank that holds the cell, carries tile-local
// indices with the tile id that makes them resolvable, and is rebuilt for free
// on every regrid because flagfield already runs there.
//
// An entry is emitted for each in-domain cell face whose neighbour is outside
// the fluid (flag4 < 0). The face code follows gcdf_update's convention:
//
//     i-1 -> 1 (X_NEG)    i+1 -> 4 (X_POS)
//     j-1 -> 3 (Y_NEG)    j+1 -> 2 (Y_POS)
//     k-1 -> 5 (Z_NEG)    k+1 -> 6 (Z_POS)
//
// which makes the face->bcside map an identity: face code cs takes bcside<cs>
// (x-/x+ = 1/4, y-/y+ = 3/2, z-/z+ = 5/6, as documented in gc_flagfield). A
// neighbour that is outside the DOMAIN gets that bcside; an interior solid/topo
// neighbour gets WALL. patchBC_gcb_convert still overrides bc afterwards, as it
// did for the file-derived list.
//
// Periodic directions need no special case: FillBoundary fills the exterior
// ghost from the wrapped interior, so the neighbour reads > 0 and no entry is
// produced.
//
// Ordering note — this must run after flag4's -1 -> OBJ_FLAG conversion and its
// FillBoundary (so the domain-exterior ghost ring reads OBJ_FLAG), and before
// flagfield's GC4LOOP tagging of flag1/2/3, which consumes what this produces.
// ghostcell::flagfield calls it at exactly that point.
// ---------------------------------------------------------------------------

void ghostcell::gcb4_generate(lexer *p)
{
    const int WALL_BC = 21;

    // Face code cs -> bcside<cs>. Index 0 unused.
    int bcs[7] = {0, p->bcside1, p->bcside2, p->bcside3,
                     p->bcside4, p->bcside5, p->bcside6};

    if(p->B98>=3)
    {
        for(auto &e : bcs)
        {
            if(e == 6)
            {
                e = 1;
            }
        }
    }

    // Domain bounds in the global index space, and the tile origin that turns a
    // tile-local index into a global one: GLOBAL_I/J/K, LEVEL_DOMAIN_DECL and
    // DOMAIN_LO/HI, from global_index.h.

    // -----------------------------------------------------------------------
    // count
    #if USE_AMREX
    const int nlevels = p->nlevs;
    #else
    const int nlevels = 1;
    #endif

    p->gcb4.resize_levels(nlevels);

    LEVEL_LOOP
    {
        count=0;
        TILE_LOOP
        IJKLOOP
        PCHECK
        {
            if(p->flag4(i-1,j,k)<0) ++count;
            if(p->flag4(i+1,j,k)<0) ++count;
            if(p->flag4(i,j-1,k)<0) ++count;
            if(p->flag4(i,j+1,k)<0) ++count;
            if(p->flag4(i,j,k-1)<0) ++count;
            if(p->flag4(i,j,k+1)<0) ++count;
        }

        const int newcount = count;

        // -----------------------------------------------------------------------
        // size the lists
        p->gcb4[p->level].assign(newcount, {});

        // -----------------------------------------------------------------------
        // fill

        LEVEL_DOMAIN_DECL(dom)

        count = 0;
        TILE_LOOP
        {
            IJKLOOP
            PCHECK
            {
                const int gi = GLOBAL_I;
                const int gj = GLOBAL_J;
                const int gk = GLOBAL_K;

                // cs, neighbour flag, and whether that neighbour is outside the domain
                const int  cs_list[6]  = {1, 4, 3, 2, 5, 6};
                const int  nb_flag[6]  = {p->flag4(i-1,j,k), p->flag4(i+1,j,k),
                                        p->flag4(i,j-1,k), p->flag4(i,j+1,k),
                                        p->flag4(i,j,k-1), p->flag4(i,j,k+1)};
                const bool nb_out[6]   = {gi-1 < DOMAIN_LO(dom,0), gi+1 > DOMAIN_HI(dom,0),
                                        gj-1 < DOMAIN_LO(dom,1), gj+1 > DOMAIN_HI(dom,1),
                                        gk-1 < DOMAIN_LO(dom,2), gk+1 > DOMAIN_HI(dom,2)};

                for(int d=0; d<6; ++d)
                {
                    if(nb_flag[d] >= 0)
                        continue;

                    auto &e = p->gcb4[p->level][count];

                    e.i  = i;
                    e.j  = j;
                    e.k  = k;
                    e.cs = cs_list[d];
                    e.bc = nb_out[d] ? bcs[cs_list[d]] : WALL_BC;

                    // Ghost-cell distance: half the cell spacing normal to the face.
                    // d/2 is the axis (0,1 -> x; 2,3 -> y; 4,5 -> z), matching cs_list.
                    //
                    // Set here rather than in a later pass because this is the only point
                    // at which the entry's tile context and level are already installed:
                    // DXP[IP] resolves through ORIGIN_I, which folds in both amr_tile_lo
                    // and max_i*level, so the spacing read is this level's. A separate
                    // pass would have to replay the context per entry to get the same
                    // answer, and would then be undone by the next regrid's regeneration.
                    e.dist = 0.5*(d<2 ? p->DXP[IP] : (d<4 ? p->DYP[JP] : p->DZP[KP]));

                    GCB_SET_TILE(e);

                    ++count;
                }
            }
        }
        assert(count == newcount);

        if(p->level == 0)
            p->gcb4_count = newcount;
    }

}
