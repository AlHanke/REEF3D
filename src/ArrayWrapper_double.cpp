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

#include "ArrayWrapper_double.h"
#include "lexer.h"

ArrayWrapper_double::ArrayWrapper_double(lexer* p) : p(p)
{
}

void ArrayWrapper_double::resize(double default_value)
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

double& ArrayWrapper_double::operator[] (int index)
{
    #if USE_AMREX
    const int jk = p->jmax * p->kmax;
    const int i = index / jk;
    const int rem = index - i * jk;
    const int j = rem / p->kmax;
    const int k = rem - j * p->kmax;

    const auto lo = amrex::lbound(p->amr_cell_mfi->tilebox());
    const int ii = lo.x + i;
    const int jj = lo.y + j;
    const int kk = lo.z + k;

    return data[p->level][*(p->amr_cell_mfi)].array()(ii, jj, kk, 0);
    #else
    return data[p->level][index];
    #endif
}

ArrayWrapper_double::operator double* ()
{
    #if USE_AMREX
    return data[p->level][*(p->amr_cell_mfi)].dataPtr(0);
    #else
    return data[p->level].data();
    #endif
}

void ArrayWrapper_double::setVal(double val, bool includeGhost)
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
void ArrayWrapper_double::fillBoundary()
{
    LEVEL_LOOP
    data[p->level].FillBoundary();
}

void ArrayWrapper_double::fillHigherLevels()
{
    const amrex::IntVect ratio = p->ref_vec;

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const auto& crse_mf = data[lev-1];
        auto& fine_mf = data[lev];

        for (amrex::MFIter fmfi(fine_mf); fmfi.isValid(); ++fmfi)
        {
            const amrex::Box& fine_valid_box = fmfi.validbox();
            const amrex::Box coarse_target_box = amrex::coarsen(fine_valid_box, ratio);

            auto const& fine_arr = fine_mf.array(fmfi);

            for (amrex::MFIter cmfi(crse_mf); cmfi.isValid(); ++cmfi)
            {
                const amrex::Box overlap_coarse_box = coarse_target_box & cmfi.validbox();
                if (!overlap_coarse_box.ok())
                {
                    continue;
                }

                const amrex::Box overlap_fine_box = amrex::refine(overlap_coarse_box, ratio) & fine_valid_box;
                auto const& crse_arr = crse_mf.const_array(cmfi);
                const int ref_ratio = p->ref_ratio;

                amrex::ParallelFor(overlap_fine_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int ic = amrex::coarsen(i, ref_ratio);
                    const int jc = amrex::coarsen(j, ref_ratio);
                    const int kc = amrex::coarsen(k, ref_ratio);
                    fine_arr(i, j, k, 0) = crse_arr(ic, jc, kc, 0);
                });
            }
        }

        fine_mf.FillBoundary();
    }
}
#endif
