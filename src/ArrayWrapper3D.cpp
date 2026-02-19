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

#include "ArrayWrapper3D.h"
#include "ArrayWrapper3D_imp.h"
#include "lexer.h"
#include <algorithm>

ArrayWrapper3D::ArrayWrapper3D(lexer *pp, unsigned int _data_location) : p(pp), data_location(_data_location)
{
}

ArrayWrapper3D::~ArrayWrapper3D()
{
    // out-of-line: the USE_AMREX path gives this a body (deregister_imf)
}

void ArrayWrapper3D::resize(int default_value)
{
    #if USE_AMREX
    data.resize(p->nlevs);
    LEVEL_LOOP
    {
        data[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], 1, p->margin);
        data[p->level].setVal(default_value);
    }
    #else
    // Single level: one flat array, sized from the same lexer metrics the IJK
    // macro reads. grid::assign_margin sets imax/jmax/kmax and imin/jmin/kmin
    // together, so both are final by the time resize runs.
    data.resize(static_cast<std::size_t>(p->imax)*p->jmax*p->kmax, default_value);
    cache_addressing();
    #endif
}

void ArrayWrapper3D::setVal(int val, bool includeGhost)
{
    #if USE_AMREX
    LEVEL_LOOP
    {
        data[p->level].setVal(val, (includeGhost ? p->margin : 0));
    }
    #else
    if(includeGhost)
    {
        std::fill(data.begin(), data.end(), val);
    }
    else
    {
        int i,j,k;
        LOOP
        {
            operator()(i,j,k) = val;
        }
    }
    #endif
}

ArrayWrapper3D::operator int *()
{
    #if USE_AMREX
    return data[p->level][*(p->amr_cell_mfi)].dataPtr(0);
    #else
    return data.data(); // unshifted: callers index this with IJK, not with i/j/k
    #endif
}

#if USE_AMREX
void ArrayWrapper3D::fillBoundary()
{
    LEVEL_LOOP
    data[p->level].FillBoundary();
}

void ArrayWrapper3D::fillHigherLevels()
{
    const amrex::IntVect ratio = p->ref_ratio * amrex::IntVect::TheUnitVector();
    const int ratio_x = ratio[0];
    const int ratio_y = ratio[1];
    const int ratio_z = ratio[2];

    int dir = -1;
    if (data_location == 1) dir = 0;
    else if (data_location == 2) dir = 1;
    else if (data_location == 3) dir = 2;

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const auto& crse_mf = data[lev-1];
        auto& fine_mf = data[lev];

        amrex::BoxArray coarsened_fine_ba = amrex::coarsen(fine_mf.boxArray(), ratio);

        amrex::iMultiFab coarse_on_fine_layout(coarsened_fine_ba, fine_mf.DistributionMap(), 1, 0);
        coarse_on_fine_layout.ParallelCopy(crse_mf, 0, 0, 1);

        for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& fine_valid_box = mfi.validbox();
            auto const& fine_arr = fine_mf.array(mfi);
            auto const& crse_arr = coarse_on_fine_layout.const_array(mfi);

            amrex::ParallelFor(fine_valid_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int ic = amrex::coarsen(i, ratio_x);
                const int jc = amrex::coarsen(j, ratio_y);
                const int kc = amrex::coarsen(k, ratio_z);
                fine_arr(i, j, k, 0) = crse_arr(ic, jc, kc, 0);
            });

            if(dir != -1)
            {
                if(fine_valid_box.bigEnd(dir) == p->amrex_geometry[lev].Domain().bigEnd(dir))
                {
                    amrex::Box end_box = fine_valid_box;
                    end_box.setBig(dir, p->amrex_geometry[lev].Domain().bigEnd(dir)-1);
                    end_box.setSmall(dir, p->amrex_geometry[lev].Domain().bigEnd(dir)-2);

                    int ii = 0, jj = 0, kk = 0;
                    if(dir==0) ii = -1;
                    else if(dir==1) jj = -1;
                    else if(dir==2) kk = -1;
                    amrex::ParallelFor(end_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        fine_arr(i, j, k, 0) = fine_arr(i+ii, j+jj, k+kk, 0);
                    });
                }
            }
        }

        fine_mf.FillBoundary();
    }
}

#else

/*!
    * @brief Precomputes the flat addressing used by operator() and operator[].
    *
    * The origin and both strides are folded into m_base, so that
    *
    *   data[(i-imin)*jmax*kmax + (j-jmin)*kmax + (k-kmin)]  ==  m_base[i*m_js + j*m_ks + k]
    *
    * and operator() needs no lexer dereference at all. That is what lets the
    * index be shared with neighbouring accesses instead of being rebuilt from
    * p on every one — the wrapper's own p is a different pointer to the
    * compiler than the p at the call site, so nothing computed through it can
    * be reused.
    *
    * m_base is not a pointer formed before the array: imin/jmin/kmin are all
    * -margin, so it lands at data.data() + margin*(m_js + m_ks + 1).
    *
    * The strides are long deliberately. An int store through the reference
    * returned by operator() may alias an int member under TBAA, which forces a
    * reload of the strides on every iteration of a writing loop; as longs they
    * stay in registers and the address strength-reduces to an induction
    * variable.
*/
void ArrayWrapper3D::cache_addressing() noexcept
{
    m_ks   = p->kmax;
    m_js   = static_cast<long>(p->jmax) * p->kmax;
    m_base = data.data() - p->imin*m_js - p->jmin*m_ks - p->kmin;
}
#endif
