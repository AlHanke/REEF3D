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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include "wsf_locate.h"
#include "lexer.h"

int wsf_locate::cell_1d(const std::vector<double>& N, int org, int imax, double s)
{
    if(imax<0)
    return -1;

    // Half-open on the right, so a coordinate sitting exactly on an interior
    // node belongs to the cell above it and is claimed by exactly one tile.
    if(s<N[org] || s>=N[org+imax+1])
    return -1;

    int lo=0, hi=imax;

    while(lo<hi)
    {
        const int mid = lo + (hi-lo+1)/2;

        if(s>=N[org+mid])
        lo = mid;

        else
        hi = mid-1;
    }

    return lo;
}

bool wsf_locate::uncovered(lexer *p, int i, int j, int k)
{
    #if USE_AMREX
    // amr_cell_mf convention (grid_amrex.cpp): the finest level is a blanket 0,
    // coarser levels are makeFineMask(crse_value=1, fine_value=0) -- so
    // uncovered==1 and covered==0. Testing mask!=0 without the level guard
    // would reject every cell of the finest level.
    if(p->level >= p->nlevs-1)
    return true;

    const auto mask = p->amr_cell_mf[p->level].const_array(p->amr_fab_mfi_idx);

    // Array4 is addressed with GLOBAL indices; inside TILE_LOOP (i,j,k) are
    // offsets from amr_tile_lo.
    return mask(i + p->amr_tile_lo.x,
                j + p->amr_tile_lo.y,
                k + p->amr_tile_lo.z) != 0;
    #else
    return true;
    #endif
}
