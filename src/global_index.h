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

#ifndef GLOBAL_INDEX_H_
#define GLOBAL_INDEX_H_

// =====================================================================
// Global index and domain extent, in the CURRENT LEVEL's index space.
//
// Two questions come up wherever a boundary is identified from a loop index:
// "what is this cell's global index?" and "is that inside the domain?". Both
// answers differ between builds, and under AMReX both differ between levels, so
// they are collected here rather than re-derived per file. gcb4_generate and
// flagfield each had their own copy; fillgcb1-3 would have been a third.
//
//   GLOBAL_I/J/K      loop index -> global index
//   LEVEL_DOMAIN_DECL declare the level's domain box
//   DOMAIN_LO/HI      its bounds, per direction
//
// ---------------------------------------------------------------------
// What is wrong with the legacy form under AMReX
// ---------------------------------------------------------------------
// The expression these replace is `i + p->origin_i` compared against
// `p->gknox-1`. Both operands are wrong once AMReX owns the decomposition:
//
//   origin_i is the rank origin of the DIVEMesh decomposition (read_grid.cpp).
//   Under AMReX the domain is chopped into a BoxArray with its own ownership,
//   and a loop index is an offset from the installed TILE, not from that origin.
//
//   gknox is the level-0 global cell count. A refined level has its own, larger,
//   index space, so gknox-1 does not name its last cell.
//
// Hence amr_tile_lo for the origin and the LEVEL's Geometry for the extent. At
// level 0 the two agree — DOMAIN_HI(dom,0) is gknox-1 — so switching a call site
// over does not change level-0 behaviour.
//
// ---------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------
//   LEVEL_DOMAIN_DECL(dom)                       // once per level, outside the
//                                                // hot loop
//   if(GLOBAL_I < DOMAIN_LO(dom,0)) ...
//
// Requirements at the point of use:
//   - a lexer* named p, and loop variables i/j/k (the increment members), which
//     is what every loop macro in looping.h already assumes;
//   - a tile context installed — inside a TILE_LOOP, or reinstated for a stored
//     entry via GCB_TILE / GCB4_TILE / GCB_APPLY_TILE. Without one the indices
//     resolve against the default whole-box origin;
//   - p->level equal to the level being addressed, since LEVEL_DOMAIN_DECL reads
//     the geometry of p->level.
//
// LEVEL_DOMAIN_DECL declares nothing without AMReX, and DOMAIN_LO/HI ignore their
// box argument there — so a call site needs no #if of its own. The argument is
// named rather than assumed so the dependency is visible at the point of use; an
// earlier version relied on a local that happened to be called `dom`.
//
// Hoist the declaration out of the loop where the loop is hot: DOMAIN_HI expands
// to a member call on the box, and gcb4_generate evaluates it six times per cell.
// =====================================================================

#if USE_AMREX

#include <AMReX_Box.H>

#define GLOBAL_I (i + p->amr_tile_lo.x)
#define GLOBAL_J (j + p->amr_tile_lo.y)
#define GLOBAL_K (k + p->amr_tile_lo.z)

#define LEVEL_DOMAIN_DECL(name) \
    const amrex::Box name = p->amrex_geometry[p->level].Domain();

#define DOMAIN_LO(dom, d) ((dom).smallEnd(d))
#define DOMAIN_HI(dom, d) ((dom).bigEnd(d))

#else

#define GLOBAL_I (i + p->origin_i)
#define GLOBAL_J (j + p->origin_j)
#define GLOBAL_K (k + p->origin_k)

// No levels and no Geometry: the extent is the single global grid, and the box
// argument below is discarded unexpanded, so call sites need not declare one.
#define LEVEL_DOMAIN_DECL(name)

#define DOMAIN_LO(dom, d) (0)
#define DOMAIN_HI(dom, d) ((d)==0 ? p->gknox-1 : ((d)==1 ? p->gknoy-1 : p->gknoz-1))

#endif

#endif
