/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef FIELDINT_AMREX_H_
#define FIELDINT_AMREX_H_

#include "definitions.h"   // DataLocation
#include "fieldint.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

/// @p location is registered as well as used here: grid_amrex::ba_for must give
/// the container the same index type on every regrid redefine, so registering as
/// cell-centred would let the first regrid silently flatten a NODE_Z field back.
static amrex::Vector<amrex::iMultiFab> make_imf(lexer* p, int ncomp,
                                                 amrex::Vector<amrex::iMultiFab>* dest,
                                                 DataLocation location = DataLocation::CELL_CENTERED)
{
    if(location == DataLocation::NODE_Z && ncomp != 1)
    {
        amrex::Abort("make_imf: NODE_Z location is only valid 1 component fieldints");
    }

    amrex::Vector<amrex::iMultiFab> result(p->nlevs);
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        result[lev].define(p->ba_for(location, lev),
                           p->amrex_distribution_mapping[lev],
                           ncomp, p->margin);
        result[lev].setVal(0, 0, ncomp, p->margin);
    }
    p->register_imf(dest, ncomp, location);
    return result;
}

class fieldint_amrex : public fieldint
{
public:
    virtual ~fieldint_amrex();

    inline int& operator()(int ii, int jj, int kk) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, kk + m_cached_oz, 0);
    }

    inline const int& operator()(int ii, int jj, int kk) const noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, kk + m_cached_oz, 0);
    }

    inline int& operator()(const amrex::IntVect& iv, int comp = 0) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(iv, comp);
    }

    inline const int& operator()(const amrex::IntVect& iv, int comp = 0) const noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(iv, comp);
    }

    void setVal(int val, bool includeGhost = false) override final;

    void FillBoundary() override final;

    inline amrex::iMultiFab& GetMultiFab() {return mf[p->level];};
    inline const amrex::iMultiFab& GetMultiFab() const {return mf[p->level];};
    inline amrex::iMultiFab& GetMultiFab(int level) {return mf[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return mf[level];};

protected:
    fieldint_amrex(lexer* p);

    /// @p location selects the storage layout. NODE_Z is the sigma-grid
    /// vertical-node layout -- one z-plane more than there are cells -- and is
    /// the only location here with a non-cell-centred index type. Unlike
    /// fields_amrex.h's field7 this carries no boundary-condition machinery:
    /// the int fields only need the internal exchange, which FillBoundary
    /// (plus the NODE_Z OverrideSync) provides.
    fieldint_amrex(lexer *p, DataLocation location);

    lexer *p = nullptr;
    DataLocation data_location = DataLocation::CELL_CENTERED;
    amrex::Vector<amrex::iMultiFab> mf = {};

private:
    /// Refreshes the Array4<int> cache when the current tile or level changes.
    /// Tile bounds are pre-computed by TILE_LOOP via set_tile_mfi.
    /// Const so the const and non-const accessors share one cache; the cache
    /// state is mutable because filling it is logically const.
    AMREX_FORCE_INLINE void refresh_cache_if_needed() const
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
        {
            // array() on a const FabArray only yields Array4<const int>, so the fab
            // is un-consted here. Writes through the cache are only reachable via
            // the non-const operators, which require a non-const fieldint_amrex.
            m_cached_arr4    = const_cast<amrex::iMultiFab&>(mf[cur_lev])
                                   .atLocalIdx(p->amr_local_fab_idx).array();
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_oz      = p->amr_tile_lo.z;
            m_cached_mfi_idx = cur_idx;
            m_cached_level   = cur_lev;
            m_cached_til_idx = cur_tile_index;
        }
        else if(cur_tile_index != m_cached_til_idx)
        {
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_oz      = p->amr_tile_lo.z;
            m_cached_til_idx = cur_tile_index;
        }
    }

    mutable amrex::Array4<int> m_cached_arr4 = {};
    mutable int m_cached_ox      = 0;
    mutable int m_cached_oy      = 0;
    mutable int m_cached_oz      = 0;
    mutable int m_cached_mfi_idx = -1;
    mutable int m_cached_level   = -1;
    mutable int m_cached_til_idx = -1;
};

#endif
#endif
