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

#ifndef GCB_SL_LIST_H_
#define GCB_SL_LIST_H_

// =====================================================================
// gcb_list_t — a ghost-cell boundary list, one std::vector of entries per
// AMR level. Templated on the entry type, so the same container serves
// multiple dimensionalities and usecases:
//
//   gcb_sl_cs_bc_list = gcb_list_t<gcb_sl_cs_bc>  2D slice lists: gcbsl1/2/4
//   gcb_list    = gcb_list_t<gcb_field_cs_bc_row> 3D field lists: gcb1/2/3/4
//
//   gcb_sl_in_list  = gcb_list_t<gcb_sl_0>  2D slice inflow list: gcbslin
//   gcb_sl_out_list = gcb_list_t<gcb_sl_cs> 2D slice outflow list: gcbslout
//
// The container never inspects the entry, so the only thing that varies is
// the payload (gcb_field carries k, gcb_sl does not). Loop macros and any
// code that only reads .i/.j/.cs/.bc work against either.
//
// Level selection is DELIBERATELY explicit: operator[] takes the level, it
// does not read p->level itself. A gcb list is a per-level container, not a
// field, so picking the list is a caller decision — patchBC_2D_fillobj, for
// one, pins p->level=0 because patch BCs are a level-0 concept, and that has
// to stay visible at the call site. An implicit p->level here would turn a
// wrong-level access into silently mis-written ghost cells with no
// diagnostic. The GCSLB1/QGCSLB1/QQGCSLB1 macros in looping2D.h carry
// p->level for the loops that do want the current level.
//
// This header must not reference lexer: it is included into lexer.h before
// the lexer class is formed. tile_ctx.h is deliberately free of that
// dependency for the same reason.
//
// ---------------------------------------------------------------------
// Tile context on every entry
// ---------------------------------------------------------------------
// Under AMReX the i/j/k on these entries are TILE-LOCAL — offsets from the
// amr_tile_lo of the tile that recorded them. They are only meaningful while
// that tile's context is installed. Entries are produced inside a TILE_LOOP
// and consumed by loops that run outside one, where the default whole-box
// context is active, so the indices would otherwise resolve against the wrong
// origin as soon as tiling yields more than one tile per box.
//
// Each entry therefore records the dense id of the tile it was written under:
//
//   producer:  e.ctx_id = p->tile_ctx_id();          // inside the TILE_LOOP
//   consumer:  GCB_TILE(e, lev);                     // before touching i/j/k
//              ...
//              GC_TILE_RESET;                        // after the loop
//
// An id, not a TileCtx: there are far fewer tiles than entries, so the shared
// table in grid_amrex costs 48 bytes per tile while each entry pays 4. The id
// is a prefix sum over (level, LocalIndex, LocalTileIndex), so producers that
// iterate in different orders still agree on it.
//
// TILE_CTX_DEFAULT means the entry was recorded with no tile installed. That is
// not an error — tile_ctx_by_id resolves it to the level's whole-box context,
// so such an entry replays exactly as it was written — but it is worth knowing
// about: see grid_amrex::tile_ctx_default_records().
//
// Restoring a context also reinstates its level, which must agree with the
// gcb_list_t bucket the entry lives in: the list's level selects the list,
// the context's level makes the FAB indices inside it resolvable.
//
// The member is AMReX-only. It cannot be written without AMReX anyway —
// tile_ctx_id() lives on grid_amrex — so producers and consumers are already
// guarded (GCB_TILE and GC_TILE_RESET compile to nothing there), and this way
// the non-AMReX build pays no per-entry cost for a concept it does not have.
// =====================================================================

#include "tile_ctx.h"

#include <cassert>
#include <vector>

// One slice inflow ghost-cell boundary entry. Replaces the old int[2] row,
// whose [0]/[1] were i/j.
struct gcb_sl
{
    int i = 0;
    int j = 0;

    #if USE_AMREX
    int ctx_id = TILE_CTX_DEFAULT;   ///< tile that i/j are relative to
    #endif
};

// One slice outflow ghost-cell boundary entry. Replaces the old int[3] row,
// whose [0]/[1] were i/j/cs.
struct gcb_sl_cs
{
    int i  = 0;
    int j  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4

    #if USE_AMREX
    int ctx_id = TILE_CTX_DEFAULT;   ///< tile that i/j are relative to
    #endif
};

// One slice ghost-cell boundary entry. Replaces the old int[4] row, whose
// [0]/[1]/[2]/[3] were i/j/cs/bc.
struct gcb_sl_cs_bc
{
    int i  = 0;
    int j  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4
    int bc = 21;   ///< boundary condition flag; 21 = default wall

    #if USE_AMREX
    int ctx_id = TILE_CTX_DEFAULT;   ///< tile that i/j are relative to
    #endif
};

// One field ghost-cell boundary entry. Replaces the old int[5] row, whose
// [0]/[1]/[2]/[3]/[4] were i/j/k/cs/bc.
struct gcb_field_cs_bc_row
{
    int i  = 0;
    int j  = 0;
    int k  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4, Z_POS=5, Z_NEG=6
    int bc = 21;   ///< boundary condition flag; 21 = default wall
    int row = -1;   ///< row in the matrix for cval, or -1 if not used

    #if USE_AMREX
    int ctx_id = TILE_CTX_DEFAULT;   ///< tile that i/j/k are relative to
    #endif
};

template<typename ENTRY>
class gcb_list_t
{
public:
    using entry_type = ENTRY;
    using level_type = std::vector<ENTRY>;

    // Unchecked in release: the assert is the bounds check, and -DNDEBUG
    // (the release target) removes it entirely. std::vector::at() would not
    // optimize away at -O3 — the compiler cannot prove lev is in range, so
    // the compare and the throw path survive.
    level_type& operator[](int lev) noexcept
    {
        assert(lev >= 0 && lev < static_cast<int>(m_lev.size())
               && "gcb_list_t: level out of range");
        return m_lev[static_cast<size_t>(lev)];
    }

    const level_type& operator[](int lev) const noexcept
    {
        assert(lev >= 0 && lev < static_cast<int>(m_lev.size())
               && "gcb_list_t: level out of range");
        return m_lev[static_cast<size_t>(lev)];
    }

    // Number of entries on a level, as int — the loop macros compare against
    // an int counter, and the dev target builds with -Wsign-conversion.
    int ssize(int lev) const noexcept
    {
        return static_cast<int>((*this)[lev].size());
    }

    void resize_levels(int nlevs)
    {
        assert(nlevs > 0 && "gcb_list_t: need at least one level");
        m_lev.resize(static_cast<size_t>(nlevs));
    }

    int nlevels() const noexcept { return static_cast<int>(m_lev.size()); }

private:
    // Born with one empty level. A list can be legitimately empty and still be
    // read: makegrid2D_basic (FNPF) seeds only mgcslice4, yet gcsl_setbcio
    // still runs GCSL1LOOP/GCSL2LOOP over gcbsl1/gcbsl2 — which correctly find
    // nothing. This is the invariant the old Iarray(gcbsl1,1,5) in read_grid
    // provided; without it those loops index level 0 of an empty container.
    std::vector<level_type> m_lev = std::vector<level_type>(1);
};

using gcb_sl_cs_bc_list = gcb_list_t<gcb_sl_cs_bc>;   ///< 2D: gcbsl1/2/4
using gcb_list    = gcb_list_t<gcb_field_cs_bc_row>;  ///< 3D: gcb1/2/3/4

using gcb_sl_list = gcb_list_t<gcb_sl>;      ///< 2D: gcbslin
using gcb_sl_cs_list = gcb_list_t<gcb_sl_cs>;  ///< 2D: gcbslout

#endif
