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

#ifndef ARRAYWRAPPER2D_IMP_H_
#define ARRAYWRAPPER2D_IMP_H_

// Inline hot-path methods of ArrayWrapper2D. Included at the bottom of
// lexer.h, after the lexer class is complete (same pattern as
// ArrayWrapper3D_imp.h) so these bodies can dereference p->... members.

#include "ArrayWrapper2D.h"
#if USE_AMREX
#include "lexer.h"
#include <cassert>
#endif

// Origin and stride are folded into m_base by cache_addressing(), so these
// touch no lexer member. Equivalent to data[(ii-imin)*jmax + (jj-jmin)].
inline int &ArrayWrapper2D::operator()(int ii, int jj) noexcept
{
    #if USE_AMREX
    refresh_cache_if_needed();
    return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, 0, 0);
    #else
    return m_base[ii*m_js + jj];
    #endif
}

inline const int &ArrayWrapper2D::operator()(int ii, int jj) const noexcept
{
    #if USE_AMREX
    refresh_const_cache_if_needed();
    return m_cached_const_arr4(ii + m_cached_const_ox, jj + m_cached_const_oy, 0, 0);
    #else
    return m_base[ii*m_js + jj];
    #endif
}

inline int &ArrayWrapper2D::operator[](int index) noexcept
{
    #if USE_AMREX
    assert(p->level == 0
           && "ArrayWrapper2D::operator[] is level-0 only under AMReX");
    refresh_cache_if_needed();

    const int ii = index / p->jmax;
    const int jj = index - ii * p->jmax;

    return m_cached_arr4(m_cached_ox + ii + p->imin,
                         m_cached_oy + jj + p->jmin, 0, 0);
    #else
    return data.data()[index];
    #endif
}

inline const int &ArrayWrapper2D::operator[](int index) const noexcept
{
    #if USE_AMREX
    assert(p->level == 0
           && "ArrayWrapper2D::operator[] is level-0 only under AMReX");
    refresh_const_cache_if_needed();

    const int ii = index / p->jmax;
    const int jj = index - ii * p->jmax;

    return m_cached_const_arr4(m_cached_const_ox + ii + p->imin,
                               m_cached_const_oy + jj + p->jmin, 0, 0);
    #else
    return data.data()[index];
    #endif
}

#if USE_AMREX

amrex::iMultiFab& ArrayWrapper2D::GetMultiFab()
{
    return GetMultiFab(p->level);
}

const amrex::iMultiFab& ArrayWrapper2D::GetMultiFab() const
{
    return GetMultiFab(p->level);
}

AMREX_FORCE_INLINE void ArrayWrapper2D::refresh_cache_if_needed() noexcept
{
    const int cur_lev  = p->level;
    const int cur_idx  = p->amr_fab_mfi_idx;
    const int cur_tile = p->amr_local_tile_idx;
    if(cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
    {
        m_cached_arr4    = GetMultiFab(cur_lev).array(*(p->amr_cell_mfi));
        m_cached_ox      = p->amr_tile_lo.x;
        m_cached_oy      = p->amr_tile_lo.y;
        m_cached_mfi_idx = cur_idx;
        m_cached_level   = cur_lev;
        m_cached_til_idx = cur_tile;
    }
    else if(cur_tile != m_cached_til_idx)
    {
        m_cached_ox      = p->amr_tile_lo.x;
        m_cached_oy      = p->amr_tile_lo.y;
        m_cached_til_idx = cur_tile;
    }
}

AMREX_FORCE_INLINE void ArrayWrapper2D::refresh_const_cache_if_needed() const noexcept
{
    const int cur_lev  = p->level;
    const int cur_idx  = p->amr_fab_mfi_idx;
    const int cur_tile = p->amr_local_tile_idx;
    if(cur_lev != m_cached_const_level || cur_idx != m_cached_const_mfi_idx)
    {
        m_cached_const_arr4    = GetMultiFab(cur_lev).const_array(*(p->amr_cell_mfi));
        m_cached_const_ox      = p->amr_tile_lo.x;
        m_cached_const_oy      = p->amr_tile_lo.y;
        m_cached_const_mfi_idx = cur_idx;
        m_cached_const_level   = cur_lev;
        m_cached_const_til_idx = cur_tile;
    }
    else if(cur_tile != m_cached_const_til_idx)
    {
        m_cached_const_ox      = p->amr_tile_lo.x;
        m_cached_const_oy      = p->amr_tile_lo.y;
        m_cached_const_til_idx = cur_tile;
    }
}

#endif
#endif
