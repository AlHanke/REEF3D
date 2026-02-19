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

#ifndef ARRAYWRAPPER3D_IMP_H_
#define ARRAYWRAPPER3D_IMP_H_

#include "ArrayWrapper3D.h"
#if USE_AMREX
#include "lexer.h"
#endif

inline int &ArrayWrapper3D::operator() (int i, int j, int k) noexcept
{
    #if USE_AMREX
    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    return data[p->level][*(p->amr_cell_mfi)].array()(lo.x + i, lo.y + j, lo.z + k, 0);
    #else
    // Origin and strides are folded into m_base by cache_addressing(), so this
    // touches no lexer member. Equivalent to data[IJK].
    return m_base[i*m_js + j*m_ks + k];
    #endif
}

inline const int &ArrayWrapper3D::operator()(int i, int j, int k) const noexcept
{
    #if USE_AMREX
    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    return data[p->level][*(p->amr_cell_mfi)].const_array()(lo.x + i, lo.y + j, lo.z + k, 0);
    #else
    // Origin and strides are folded into m_base by cache_addressing(), so this
    // touches no lexer member. Equivalent to data[IJK].
    return m_base[i*m_js + j*m_ks + k];
    #endif
}

inline int &ArrayWrapper3D::operator[] (int index) noexcept
{
    #if USE_AMREX
    // kz() is p->kmaxF for the vertical-node layout and p->kmax otherwise,
    // matching the stride the caller's flat index was built with (FIJK vs IJK).
    const int k_max = kz();
    const int jk = p->jmax * k_max;

    const int ii_encoded = index / jk;
    const int rem = index - ii_encoded * jk;
    const int jj_encoded = rem / k_max;
    const int kk_encoded = rem - jj_encoded * k_max;

    const int i = ii_encoded + p->imin;
    const int j = jj_encoded + p->jmin;
    const int k = kk_encoded + p->kmin;

    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    const int ii = lo.x + i;
    const int jj = lo.y + j;
    const int kk = lo.z + k;

    return data[p->level][*(p->amr_cell_mfi)].array()(ii, jj, kk, 0);
    #else
    return data.data()[index];
    #endif
}

inline const int &ArrayWrapper3D::operator[] (int index) const noexcept
{
    #if USE_AMREX
    // kz() is p->kmaxF for the vertical-node layout and p->kmax otherwise,
    // matching the stride the caller's flat index was built with (FIJK vs IJK).
    const int k_max = kz();
    const int jk = p->jmax * k_max;

    const int ii_encoded = index / jk;
    const int rem = index - ii_encoded * jk;
    const int jj_encoded = rem / k_max;
    const int kk_encoded = rem - jj_encoded * k_max;

    const int i = ii_encoded + p->imin;
    const int j = jj_encoded + p->jmin;
    const int k = kk_encoded + p->kmin;

    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    const int ii = lo.x + i;
    const int jj = lo.y + j;
    const int kk = lo.z + k;

    return data[p->level][*(p->amr_cell_mfi)].const_array()(ii, jj, kk, 0);
    #else
    return data.data()[index];
    #endif
}

#endif
