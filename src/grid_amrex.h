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

#if USE_AMREX
#ifndef GRID_AMREX_H_
#define GRID_AMREX_H_

#include "definitions.h"
#include "grid.h"
#include "tile_ctx.h"

// AMReX includes
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

#include <memory>
#include <vector>
#include <utility>
#include <algorithm>

class lexer;
class ghostcell;
class fdm;
class slice_amrex_prim;

namespace amrex {
    class Geometry;
    class BoxArray;
    class DistributionMapping;
}

/*!
    * @brief grid_amrex is a class that incorporates AMReX data structures and functionalities for handling the computational grid in REEF3D.
    * It inherits from the grid class and extends it with AMReX-specific members and methods.
    * The class includes methods for defining inflow and outflow boundary areas, setting up AMReX geometry, and managing AMReX MultiFabs and iMultiFabs for cell-centered data and flags.
    * It also includes members for storing AMReX Geometry, BoxArray, and DistributionMapping for each level of the AMR hierarchy.
    * The class is designed to facilitate the integration of AMReX's capabilities for handling complex grid structures and parallel computations within the REEF3D framework.
*/
class grid_amrex : public grid
{
public:
    grid_amrex() = default;
    virtual ~grid_amrex() = default;

    void regrid_amrex_box_array_and_distribution_mapping(lexer* p, fdm* a);

    void update_cell_coordinates();

    // Registry for all owning MultiFab/iMultiFab vectors so they can be
    // resized collectively during AMR regrid.
    struct MFEntry  { amrex::Vector<amrex::MultiFab>*  mf; int ncomp; DataLocation location; };
    struct IMFEntry { amrex::Vector<amrex::iMultiFab>* mf; int ncomp; DataLocation location; };
    std::vector<MFEntry>  mf_registry;
    std::vector<IMFEntry> imf_registry;

    /// Index type implied by a registry `location` (DataLocation
    /// cast to int). Only NODE_Z (7) is non-cell: z-nodal, one plane more than
    /// there are cells. The FACE_* locations keep cell-centred storage and are
    /// re-staggered at interpolation time, so they are cell-typed here.
    static amrex::IntVect location_index_type(DataLocation location) noexcept
    {
        return location == DataLocation::NODE_Z
             ? amrex::IntVect(AMREX_D_DECL(0,0,1))
             : amrex::IntVect::TheZeroVector();
    }

    /// BoxArray a registered container with @p location must be defined on at
    /// @p lev. amrex::convert shares the underlying box-list reference and
    /// preserves box count and ordering, so a NODE_Z container keeps the same
    /// DistributionMap and the same MFIter LocalIndex as every cell-centred one
    /// — TileCtx and TILE_LOOP need no knowledge of it. Converting by the zero
    /// vector is a no-op that returns the cell-centred array unchanged.
    amrex::BoxArray ba_for(DataLocation location, int lev) const
    {
        return amrex::convert(amrex_box_array[lev], location_index_type(location));
    }

    void register_mf (amrex::Vector<amrex::MultiFab>*  mf, int ncomp, DataLocation location=DataLocation::CELL_CENTERED)
    {
        mf_registry .push_back({mf, ncomp, location});
    }
    void register_imf(amrex::Vector<amrex::iMultiFab>* mf, int ncomp, DataLocation location=DataLocation::CELL_CENTERED)
    {
        imf_registry.push_back({mf, ncomp, location});
    }

    void deregister_mf(amrex::Vector<amrex::MultiFab>* mf)
    {
        auto it = std::remove_if(mf_registry.begin(), mf_registry.end(),
                                 [mf](const MFEntry& entry) { return entry.mf == mf; });
        mf_registry.erase(it, mf_registry.end());
    }
    void deregister_imf(amrex::Vector<amrex::iMultiFab>* mf)
    {
        auto it = std::remove_if(imf_registry.begin(), imf_registry.end(),
                                 [mf](const IMFEntry& entry) { return entry.mf == mf; });
        imf_registry.erase(it, imf_registry.end());
    }

    std::vector<slice_amrex_prim*> slice_registry;
    void register_slice(slice_amrex_prim* s) { slice_registry.push_back(s); }
    void deregister_slice(slice_amrex_prim* s)
    {
        auto it = std::remove(slice_registry.begin(), slice_registry.end(), s);
        slice_registry.erase(it, slice_registry.end());
    }

    // Registry for field_amrex objects so they can extend BCRecs/m_alias on regrid
    std::vector<class field_amrex*> field_registry;
    void register_field(class field_amrex* f) { field_registry.push_back(f); }
    void extend_registered_fields(int new_nlevs);
    void deregister_field(class field_amrex* f)
    {
        auto it = std::remove(field_registry.begin(), field_registry.end(), f);
        field_registry.erase(it, field_registry.end());
    }

    // Registry for WENO non-uniform-grid objects so their precomputed tables
    // are reallocated and recalculated alongside coordinate arrays on regrid.
    std::vector<class weno3_nug_func*> weno3_registry;
    std::vector<class weno_nug_func*>  weno5_registry;

    void register_weno3(class weno3_nug_func* w) { weno3_registry.push_back(w); }
    void register_weno5(class weno_nug_func*  w) { weno5_registry.push_back(w); }
    void deregister_weno3(class weno3_nug_func* w)
    {
        auto it = std::remove(weno3_registry.begin(), weno3_registry.end(), w);
        weno3_registry.erase(it, weno3_registry.end());
    }
    void deregister_weno5(class weno_nug_func* w)
    {
        auto it = std::remove(weno5_registry.begin(), weno5_registry.end(), w);
        weno5_registry.erase(it, weno5_registry.end());
    }
    void update_registered_weno(int new_nlevs);

    void resize_registered_mf(int old_nlevs, int new_nlevs);
    void fill_registered_mf_level(int lev);
    /// Redefine all registered MultiFabs at @p lev with the current BoxArray and
    /// DistributionMapping, filling data via coarse-level interpolation.
    /// Called for existing non-zero levels whose box arrays change during regrid.
    void redefine_registered_mf_level(int lev);
    /// Rebuild view-mode aliases at @p lev for all registered field_amrex objects.
    void rebuild_registered_field_aliases_level(int lev);
    void update_cell_spacing();

    // AMReX Data structures
    amrex::Vector<amrex::Geometry> amrex_geometry; // Phyiscal domain and coordinate system
    amrex::Vector<amrex::BoxArray> amrex_box_array; // BoxArray defines the index space decomposition of the domain into boxes
    amrex::Vector<amrex::DistributionMapping> amrex_distribution_mapping; // DistributionMapping defines the mapping of boxes in the BoxArray to MPI ranks
    amrex::Vector<amrex::BoxArray> amrex_slice_box_array; // BoxArray for slice planes
    amrex::Vector<amrex::DistributionMapping> amrex_slice_distribution_mapping; // DistributionMapping for slice planes
    amrex::Vector<amrex::Vector<std::pair<amrex::RealVect,amrex::RealVect>>> amrex_refined_grid_coords; // Input: Coordinates of the refined grid boxes for each level, index is offset by 1 (i.e. amrex_refined_grid_coords[0] is for level 1, etc.)

    // -----------------------------------------------------------------
    // Slice column ownership — shared by every ArrayWrapper2D.
    //
    // The int slices flatten the 3D BoxArray onto z-slabs, so every 3D box in a
    // column lands on the same (i,j) and the resulting VIEW layout overlaps.
    // Reductions over it — and any "visit each column once" loop — need a single
    // designated owner box per column.
    //
    // For an int flag the payload is column-constant, so the owner can be purely
    // geometric: the box touching the domain z-lo face. That makes the mask a
    // function of (BoxArray, DistributionMapping, Domain) alone — identical for
    // flagslice1/2/4 and IOSL — so it is built here once per grid generation
    // instead of once per instance per regrid.
    //
    // NOT shared with slice_amrex's mask: that one is ksurf-based and follows the
    // free surface, so it is data-dependent and has to be rebuilt whenever the
    // surface moves. Same layout, different predicate — they must not be merged.
    //
    // slice_view_ba is the view layout itself (one flattened box per 3D box, on
    // the 3D DistributionMapping), published so the slices define their m_view
    // against the same BoxArray the mask uses — an MFIter over m_view then
    // indexes slice_owner_mf directly.
    amrex::Vector<amrex::BoxArray>  slice_view_ba;
    amrex::Vector<amrex::iMultiFab> slice_owner_mf;

    /// Rebuild both, for all levels. Must run after the box arrays are in place
    /// and before anything that reads them — at setup, and on regrid ahead of the
    /// slice_registry callbacks, which reduce through the mask.
    void build_slice_owner();

    const amrex::BoxArray& slice_view_boxarray(int lev) const
    {
        AMREX_ASSERT(lev >= 0 && lev < int(slice_view_ba.size()));
        return slice_view_ba[lev];
    }

    /// 1 on the one view cell per column that counts, 0 on the duplicates.
    const amrex::iMultiFab& slice_owner(int lev) const
    {
        AMREX_ASSERT(lev >= 0 && lev < int(slice_owner_mf.size()));
        AMREX_ASSERT(slice_owner_mf[lev].ok() && "build_slice_owner not run for this grid");
        return slice_owner_mf[lev];
    }

    // Looping structures
    amrex::Vector<amrex::iMultiFab> amr_cell_mf;
    std::unique_ptr<amrex::MFIter> default_cell_mfi;

    // Tile bounds — set once per tile by set_tile_mfi; shared with field_amrex cache.
    amrex::Dim3 amr_tile_lo    = {0, 0, 0};
    amrex::Dim3 amr_tile_hi    = {0, 0, 0};
    int amr_fab_mfi_idx        = -1; // Index of current MFIter
    int amr_local_fab_idx      = -1; // Local FAB index within the current MFIter
    int amr_local_tile_idx     = -1; // Local tile index within the current MFIter
    int amr_tile_ctx_id        = TILE_CTX_DEFAULT; // Dense id of the installed tile

    /// Whether the installed tile owns its (i,j) columns — see TileCtx::slice_owner
    /// for the predicate. This is what PSLICEOWNER / SLICEOWNLOOP* test, so that a
    /// slice loop emitting one record per column visits each column exactly once
    /// even when a box is split into several tiles in z.
    ///
    /// Only meaningful inside a TILE_LOOP or under a context restored by
    /// set_tile_ctx. Outside one the default context is installed, which is box 0's
    /// — so on a rank owning several boxes the bit then answers for box 0 alone.
    bool amr_slice_owner       = false;

    // Tile context table, rebuilt by build_tile_ctx_table on every decomposition
    // change. tile_ctx_table is indexed by the dense id; the offsets are the
    // prefix sums that turn (level, LocalIndex, LocalTileIndex) into that id.
    std::vector<TileCtx>           tile_ctx_table;
    std::vector<TileCtx>           tile_ctx_default;      // per level, untiled whole box
    std::vector<int>               tile_ctx_level_offset; // size nlevs+1
    std::vector<std::vector<int>>  tile_ctx_fab_offset;   // [lev][local fab] -> first id

    mutable long amr_default_ctx_records = 0; // see tile_ctx_default_records()

    /// Snapshot of a tile's addressing context — defined in tile_ctx.h, which is
    /// shared with the containers that carry contexts (see gcb_sl_list.h) and so
    /// cannot depend on this header. Aliased here because the capture/restore API
    /// lives on this class and callers name it as lexer::TileCtx.
    using TileCtx = ::TileCtx;

    /// Bumped on every regrid; stamped into each TileCtx so that snapshots which
    /// outlived the decomposition they were taken from can be detected.
    int amr_grid_gen = 0;

    /// Dense id of the currently installed tile, or TILE_CTX_DEFAULT when none
    /// is (i.e. outside any TILE_LOOP). This is the value a producer stores per
    /// entry; see tile_ctx_by_id for the consumer side.
    ///
    /// Read, never derived. An id computed from amr_local_fab_idx /
    /// amr_local_tile_idx would silently alias the untiled default context onto
    /// tile 0 of box 0 — see TILE_CTX_DEFAULT in tile_ctx.h.
    inline int tile_ctx_id() const noexcept
    {
        if (amr_tile_ctx_id == TILE_CTX_DEFAULT)
            ++amr_default_ctx_records;
        return amr_tile_ctx_id;
    }

    /// Resolve a stored id. TILE_CTX_DEFAULT yields the level's untiled
    /// whole-box context, so consumers need no special case for entries that
    /// were recorded outside a tile loop.
    inline const TileCtx& tile_ctx_by_id(int id, int lev) const noexcept
    {
        if (id == TILE_CTX_DEFAULT)
        {
            AMREX_ASSERT(lev >= 0 && lev < int(tile_ctx_default.size()));
            return tile_ctx_default[size_t(lev)];
        }

        AMREX_ASSERT(id >= 0 && id < int(tile_ctx_table.size()));
        return tile_ctx_table[size_t(id)];
    }

    /// How many times tile_ctx_id() has handed out TILE_CTX_DEFAULT — i.e. how
    /// often a producer recorded indices while no tile context was installed.
    /// Those indices are relative to box 0 of the rank and are only addressable
    /// at all if the rank owns one box; replaying the context reproduces the
    /// producer's addressing faithfully, which is self-consistent but not
    /// correct. A non-zero count is a list of sites to audit, not a failure.
    long tile_ctx_default_records() const noexcept { return amr_default_ctx_records; }
    void reset_tile_ctx_default_records() noexcept { amr_default_ctx_records = 0; }

    /// Capture the currently active tile context for later restore.
    inline TileCtx tile_ctx() const noexcept
    {
        TileCtx c;
        c.lo             = amr_tile_lo;
        c.hi             = amr_tile_hi;
        c.level          = level;
        c.local_fab_idx  = amr_local_fab_idx;
        c.fab_idx        = amr_fab_mfi_idx;
        c.local_tile_idx = amr_local_tile_idx;
        c.gen            = amr_grid_gen;
        c.id             = amr_tile_ctx_id;
        c.slice_owner    = amr_slice_owner;
        return c;
    }

    /// Extract the tile context an MFIter currently points at. The level is the
    /// ambient one: an MFIter carries no level of its own, TILE_LOOP runs it
    /// inside the level set by LEVEL_LOOP.
    ///
    /// The id is a prefix sum over this rank's FABs, so it is a pure function of
    /// (level, LocalIndex, LocalTileIndex) — two producers iterating in
    /// different orders still agree. Valid only for an MFIter tiled the same way
    /// the table was built (TilingIfNotGPU); an untiled iterator gets
    /// TILE_CTX_DEFAULT, which is what default_cell_mfi must resolve to.
    inline TileCtx tile_ctx(const amrex::MFIter& mfi) const noexcept
    {
        const amrex::Box tb = mfi.tilebox();

        TileCtx c;
        c.lo             = amrex::lbound(tb);
        c.hi             = amrex::ubound(tb);
        c.level          = level;
        c.local_fab_idx  = mfi.LocalIndex();
        c.fab_idx        = mfi.index();
        c.local_tile_idx = mfi.LocalTileIndex();
        c.gen            = amr_grid_gen;
        c.id             = tile_ctx_id_of(mfi);
        c.slice_owner    = slice_owner_of(mfi, tb);
        return c;
    }

    /// The TileCtx::slice_owner predicate for an MFIter's current tile. Split out
    /// of tile_ctx(mfi) only to keep the two conditions and their guard readable.
    ///
    /// mfi.index() is used as a subscript into amrex_box_array[level]: valid for
    /// amr_cell_mf, for every registered field (they share the level's BoxArray)
    /// and for slice_view_ba (one flattened box per 3D box, in order). Anything
    /// else — or a call before the box arrays exist, which happens while the grid
    /// is still being set up — yields false rather than a guess.
    inline bool slice_owner_of(const amrex::MFIter& mfi,
                               const amrex::Box& tb) const noexcept
    {
        if (level < 0 || level >= int(amrex_box_array.size())
                      || level >= int(amrex_geometry.size()))
            return false;

        const amrex::BoxArray& ba3d = amrex_box_array[size_t(level)];

        if (mfi.index() < 0 || mfi.index() >= int(ba3d.size()))
            return false;

        const amrex::Box b3  = ba3d[mfi.index()];
        const int        zlo = amrex_geometry[size_t(level)].Domain().smallEnd(2);

        return (b3.smallEnd(2) == zlo)                    // (a) owner box
            && (amrex::lbound(tb).z == b3.smallEnd(2));   // (b) z-lowest tile of it
    }

    /// Prefix-sum id of an MFIter's current tile, or TILE_CTX_DEFAULT if the
    /// table is not built yet (early setup) or the iterator is untiled.
    inline int tile_ctx_id_of(const amrex::MFIter& mfi) const noexcept
    {
        if (level < 0 || level >= int(tile_ctx_fab_offset.size()))
            return TILE_CTX_DEFAULT;

        const auto& fab_off = tile_ctx_fab_offset[size_t(level)];
        const int   li      = mfi.LocalIndex();

        if (li < 0 || li + 1 >= int(fab_off.size()))
            return TILE_CTX_DEFAULT;

        // An untiled iterator reports 1 tile for a FAB the table split into
        // several; it is not a member of the enumeration and must not get an id.
        if (mfi.numLocalTiles() != fab_off[size_t(li) + 1] - fab_off[size_t(li)])
            return TILE_CTX_DEFAULT;

        return tile_ctx_level_offset[size_t(level)] + fab_off[size_t(li)]
             + mfi.LocalTileIndex();
    }

    /// Enumerate every tile of every level and assign dense ids, plus one
    /// untiled whole-box default context per level. Must be rebuilt whenever the
    /// decomposition changes — call sites are the initial grid setup and the end
    /// of regrid_amrex_box_array_and_distribution_mapping, alongside the
    /// amr_grid_gen bump.
    void build_tile_ctx_table();

    /// The single writer of the tile addressing state. Every field accessor
    /// resolves (i,j,k) through these members, so a context established by any
    /// route must set all of them — keeping one writer is what stops a later
    /// added member from being set on one path and left stale on the other.
    ///
    /// Deliberately does NOT write 'level': level is loop state owned by
    /// LEVEL_LOOP, not part of a tile's addressing, and TILE_LOOP's guard must be
    /// able to restore the tile geometry on exit without disturbing it. Only
    /// set_tile_ctx — reinstating a context detached from its loop — sets it.
    inline void apply_tile_ctx(const TileCtx& c) noexcept
    {
        amr_tile_lo        = c.lo;
        amr_tile_hi        = c.hi;
        amr_local_fab_idx  = c.local_fab_idx;
        amr_fab_mfi_idx    = c.fab_idx;
        amr_local_tile_idx = c.local_tile_idx;
        amr_tile_ctx_id    = c.id;
        amr_slice_owner    = c.slice_owner;
    }

    /// Install the tile an MFIter currently points at, returning the context it
    /// displaced so a guard can restore it. Called exactly once per tile loop
    /// iteration by TILE_LOOP. The returned context is a value, so it stays valid
    /// after the MFIter that produced it has advanced or gone out of scope.
    inline TileCtx set_tile_mfi(amrex::MFIter* mfi) noexcept
    {
        const TileCtx old = tile_ctx();
        apply_tile_ctx(tile_ctx(*mfi));
        return old;
    }

    /// Re-establish a captured tile context outside of TILE_LOOP, so that
    /// tile-local (i,j,k) recorded under that tile resolve correctly again.
    /// Sets 'level' too: local_fab_idx indexes a per-level FAB vector, so a
    /// context applied at the wrong level silently reads a different level's FAB.
    ///
    /// Unlike set_tile_mfi, this reinstates a context that may be arbitrarily old,
    /// so it is the path that has to check the context is still addressable.
    inline void set_tile_ctx(const TileCtx& c, const int lev) noexcept
    {
        AMREX_ASSERT(c.gen == amr_grid_gen);
        AMREX_ASSERT(c.level == lev);
        AMREX_ASSERT(c.local_fab_idx >= 0);

        apply_tile_ctx(c);
    }

    /// Restore the default whole-box context after a set_tile_ctx sequence.
    inline void reset_tile_ctx() noexcept
    {
        level = 0;
        set_tile_mfi(default_cell_mfi.get());
    }

    int nlevs = 1; // Number of AMR levels
    amrex::IntVect ref_vec;
    static constexpr int ref_ratio = 2;
    const int ncomp = 1;
    int bc_type[6] = {0,0,0,0,0,0};

    bool changed = false;

protected:
    void setup_amrex_geometry(lexer* p, ghostcell* pgc);

    static constexpr int max_nlevs = 2;
private:
    void create_amrex_box_array_and_distribution_mapping_level_n();
    void create_slice_BoxArray_and_DistributionMapping(int lev);
    void output_amrex_level_info();
};

#endif
#endif
