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

#ifndef TILE_CTX_H_
#define TILE_CTX_H_

// =====================================================================
// TileCtx — a snapshot of one tile's addressing context.
//
// Under AMReX every (i,j,k) produced inside a TILE_LOOP is TILE-LOCAL: it is
// an offset from that tile's amr_tile_lo. Code that records such indices and
// dereferences them later runs outside the loop, where a different tile (or
// the whole-box default) is installed, and the indices then resolve against
// the wrong origin. Storing a TileCtx alongside the indices and reinstating it
// via lexer::set_tile_ctx before the dereference is what makes them mean again
// what they meant when they were written.
//
// An amrex::MFIter cannot serve this purpose: it is non-copyable, owns an
// iteration position, and dies with its loop. Everything the field accessors
// actually need is the POD captured here.
//
// This header deliberately depends on nothing but AMReX_Dim3.H. It is included
// both by grid_amrex.h (which produces contexts) and by containers such as
// gcb_sl_list.h that merely carry them, and the latter are pulled into lexer.h
// before the lexer class is formed — so anything heavier would be a cycle.
//
// See grid_amrex.h for the capture/restore API: tile_ctx(), set_tile_mfi(),
// apply_tile_ctx(), set_tile_ctx().
// =====================================================================

/// Reserved id meaning "no tile context installed" — the untiled whole-box
/// default that is active outside any TILE_LOOP. It is NOT a member of the
/// tiled enumeration the id table indexes: default_cell_mfi is built with
/// tiling off, so its LocalTileIndex() is 0 and its bounds are the whole valid
/// box. Deriving an id arithmetically outside a tile loop would therefore
/// produce a valid-looking id naming tile 0 of box 0, whose lo is not the
/// box's lo. Hence ids are maintained as state, and this sentinel is what
/// producers record when they ran outside a tile loop.
constexpr int TILE_CTX_DEFAULT = -1;

#if USE_AMREX

#include <AMReX_Dim3.H>

struct TileCtx
{
    /// Tile bounds; lo is the origin that tile-local (i,j,k) are offsets from.
    amrex::Dim3 lo = {0, 0, 0};
    amrex::Dim3 hi = {0, 0, 0};

    /// AMR level. Not redundant with a per-level container's own level index:
    /// local_fab_idx below indexes a per-level FAB vector, so restoring a
    /// context at the wrong level silently reads a different level's FAB. A
    /// context stored in a per-level list must agree with the list's level.
    int level = 0;

    /// MFIter::LocalIndex() — slot in the owning rank's FAB vector. This is the
    /// handle the Array4 fetch uses. Valid only (a) on the rank that produced
    /// it, (b) at 'level', and (c) for FabArrays sharing
    /// amrex_distribution_mapping[level]. (c) holds for every FabArray
    /// reachable through the field accessors. Never MPI-communicate a TileCtx —
    /// fab_idx is the rank-independent identity if that is ever needed.
    int local_fab_idx = -1;

    /// MFIter::index() — global BoxArray subscript of the FAB. Used as the
    /// Array4 cache key, and as a stable identity in asserts and diagnostics.
    int fab_idx = -1;

    /// MFIter::LocalTileIndex() — which tile within the FAB. Under tiling
    /// neither fab_idx nor local_fab_idx is unique per iteration; the pair
    /// (fab_idx, local_tile_idx) is.
    int local_tile_idx = -1;

    /// Grid generation, see grid_amrex::amr_grid_gen. A context captured before
    /// a regrid addresses a decomposition that no longer exists; set_tile_ctx
    /// asserts on a mismatch rather than reading a stale box.
    int gen = -1;

    /// Dense index into grid_amrex's tile context table, or TILE_CTX_DEFAULT.
    /// This is what containers should store per entry (4 bytes) instead of a
    /// whole TileCtx (48) — there are far fewer tiles than entries. Carried on
    /// the struct itself so that apply_tile_ctx, the single writer of the
    /// addressing state, keeps the installed id consistent with the installed
    /// bounds; deriving one from the other after the fact is what goes wrong.
    int id = TILE_CTX_DEFAULT;
};

#else

// Without AMReX there is no tiling and (i,j,k) are already rank-global, so a
// context carries no information. Kept as a type so that containers can hold a
// member unconditionally; mark such members [[no_unique_address]] to keep them
// zero-cost in this build.
struct TileCtx {};

#endif

#endif
