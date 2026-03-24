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

#include "ArrayWrapper_int.h"
#include "lexer.h"

ArrayWrapper_int::ArrayWrapper_int(lexer* p) : p(p)
{
}

void ArrayWrapper_int::resize(int default_value)
{
    #if USE_AMREX
    data.resize(p->nlevs);
    #else
    data.resize(1);
    #endif
    LEVEL_LOOP
    {
        #if USE_AMREX
        data[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], 1, p->margin);
        data[p->level].setVal(default_value);
        #else
        // Initialize the vector with the appropriate size based on lexer parameters
        const int size = p->imax*p->jmax*p->kmax;
        data[p->level].resize(size,default_value);
        #endif
    }
};

int& ArrayWrapper_int::operator[] (int index)
{
    #if USE_AMREX
    const int jk = p->jmax * p->kmax;

    const int ii_encoded = index / jk;
    const int rem = index - ii_encoded * jk;
    const int jj_encoded = rem / p->kmax;
    const int kk_encoded = rem - jj_encoded * p->kmax;

    const int i = ii_encoded + p->imin;
    const int j = jj_encoded + p->jmin;
    const int k = kk_encoded + p->kmin;

    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    const int ii = lo.x + i;
    const int jj = lo.y + j;
    const int kk = lo.z + k;

    return data[p->level][*(p->amr_cell_mfi)].array()(ii, jj, kk, 0);
    #else
    return data[p->level][index];
    #endif
}

ArrayWrapper_int::operator int* ()
{
    #if USE_AMREX
    return data[p->level][*(p->amr_cell_mfi)].dataPtr(0);
    #else
    return data[p->level].data();
    #endif
}

void ArrayWrapper_int::setVal(int val, bool includeGhost)
{
    #if USE_AMREX
    LEVEL_LOOP
    {
        data[p->level].setVal(val, (includeGhost ? p->margin : 0));
    }
    #else
    p->level = 0;
    std::fill(data[p->level].begin(), data[p->level].end(), val);
    #endif
}

#if USE_AMREX
void ArrayWrapper_int::fillBoundary()
{
    LEVEL_LOOP
    data[p->level].FillBoundary();
}

amrex::iMultiFab& ArrayWrapper_int::GetMultiFab()
{
    return data[p->level];
}

const amrex::iMultiFab& ArrayWrapper_int::GetMultiFab() const
{
    return data[p->level];
}

amrex::iMultiFab& ArrayWrapper_int::GetMultiFab(int level)
{
    return data[level];
}

const amrex::iMultiFab& ArrayWrapper_int::GetMultiFab(int level) const
{
    return data[level];
}

void ArrayWrapper_int::fillHigherLevels()
{
    const amrex::IntVect ratio = p->ref_vec;
    const int ratio_x = ratio[0];
    const int ratio_y = ratio[1];
    const int ratio_z = ratio[2];

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
        }

        fine_mf.FillBoundary();
    }
}
#endif
