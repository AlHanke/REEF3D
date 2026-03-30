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

#include "fieldint.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

class fieldint_amrex : public fieldint
{
public:
    virtual ~fieldint_amrex() = default;

    inline int& operator()(int ii, int jj, int kk) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, kk + m_cached_oz, 0);
    }

    inline const int& operator()(int ii, int jj, int kk) const noexcept override final
    {
        refresh_const_cache_if_needed();
        return m_cached_const_arr4(ii + m_cached_const_ox, jj + m_cached_const_oy, kk + m_cached_const_oz, 0);
    }

    inline int& operator()(const amrex::IntVect& iv, int comp = 0) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(iv, comp);
    }

    inline const int& operator()(const amrex::IntVect& iv, int comp = 0) const noexcept override final
    {
        refresh_const_cache_if_needed();
        return m_cached_const_arr4(iv, comp);
    }

    void setVal(int val, bool includeGhost = false) override final;

    void FillBoundary() override final;

    inline amrex::iMultiFab& GetMultiFab() {return mf[p->level];};
    inline const amrex::iMultiFab& GetMultiFab() const {return mf[p->level];};
    inline amrex::iMultiFab& GetMultiFab(int level) {return mf[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return mf[p->level];};

protected:
    fieldint_amrex(lexer* p);

    lexer *p = nullptr;
    amrex::Vector<amrex::iMultiFab> mf = {};

private:
    /// Refreshes the Array4<int> cache when the current tile or level changes.
    /// Tile bounds are pre-computed by TILE_LOOP via set_tile_mfi.
    AMREX_FORCE_INLINE void refresh_cache_if_needed()
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
        {
            m_cached_arr4    = mf[cur_lev][*(p->amr_cell_mfi)].array();
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

    /// Refreshes the Array4<const int> cache when the current tile or level changes.
    AMREX_FORCE_INLINE void refresh_const_cache_if_needed() const
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_const_level || cur_idx != m_cached_const_mfi_idx)
        {
            m_cached_const_arr4    = mf[cur_lev][*(p->amr_cell_mfi)].const_array();
            m_cached_const_ox      = p->amr_tile_lo.x;
            m_cached_const_oy      = p->amr_tile_lo.y;
            m_cached_const_oz      = p->amr_tile_lo.z;
            m_cached_const_mfi_idx = cur_idx;
            m_cached_const_level   = cur_lev;
            m_cached_const_til_idx = cur_tile_index;
        }
        else if(cur_tile_index != m_cached_const_til_idx)
        {
            m_cached_const_ox      = p->amr_tile_lo.x;
            m_cached_const_oy      = p->amr_tile_lo.y;
            m_cached_const_oz      = p->amr_tile_lo.z;
            m_cached_const_til_idx = cur_tile_index;
        }
    }

    amrex::Array4<int> m_cached_arr4 = {};
    int m_cached_ox      = 0;
    int m_cached_oy      = 0;
    int m_cached_oz      = 0;
    int m_cached_mfi_idx = -1;
    int m_cached_level   = -1;
    int m_cached_til_idx = -1;

    mutable amrex::Array4<const int> m_cached_const_arr4 = {};
    mutable int m_cached_const_ox      = 0;
    mutable int m_cached_const_oy      = 0;
    mutable int m_cached_const_oz      = 0;
    mutable int m_cached_const_mfi_idx = -1;
    mutable int m_cached_const_level   = -1;
    mutable int m_cached_const_til_idx = -1;
};

#endif
#endif
