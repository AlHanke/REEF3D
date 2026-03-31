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

#ifndef ARRAYWRAPPER_INT_IMP_H_
#define ARRAYWRAPPER_INT_IMP_H_

#include "ArrayWrapper_int.h"
#include "lexer.h"

inline int& ArrayWrapper_int::operator[] (int index) noexcept
{
        #if USE_AMREX
        refresh_cache_if_needed();

        const int jk_max = p->jmax * p->kmax;
        const int ii_encoded = index / jk_max;
        const int rem = index - ii_encoded * jk_max;
        const int jj_encoded = rem / p->kmax;
        const int kk_encoded = rem - jj_encoded * p->kmax;

        return m_cached_arr4(m_cached_ox + ii_encoded + p->imin,
                             m_cached_oy + jj_encoded + p->jmin,
                             m_cached_oz + kk_encoded + p->kmin, 0);
        #else
        return data[p->level][index];
        #endif
    };

#if USE_AMREX

amrex::iMultiFab& ArrayWrapper_int::GetMultiFab()
{
    return data[p->level];
}

const amrex::iMultiFab& ArrayWrapper_int::GetMultiFab() const
{
    return data[p->level];
}

AMREX_FORCE_INLINE void ArrayWrapper_int::refresh_cache_if_needed() noexcept
{
    const int cur_lev = p->level;
    const int cur_idx = p->amr_fab_mfi_idx;
    const int cur_tile_index = p->amr_local_tile_idx;
    if (cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
    {
        m_cached_arr4    = data[cur_lev].array(*(p->amr_cell_mfi));
        m_cached_ox      = p->amr_tile_lo.x;
        m_cached_oy      = p->amr_tile_lo.y;
        m_cached_oz      = p->amr_tile_lo.z;
        m_cached_mfi_idx = cur_idx;
        m_cached_level   = cur_lev;
        m_cached_til_idx = cur_tile_index;
    }
    if (cur_tile_index != m_cached_til_idx)
    {
        m_cached_ox      = p->amr_tile_lo.x;
        m_cached_oy      = p->amr_tile_lo.y;
        m_cached_oz      = p->amr_tile_lo.z;
        m_cached_til_idx = cur_tile_index;
    }
}

#endif

#endif