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

#include"mgcslice2.h"
#include"lexer.h"
#include<vector>

// SLICEOWNLOOP2 rather than SLICELOOP2 — see the note in gcbslice1.cpp.
void mgcslice2::gcb_seed(lexer *p)
{
    // count gcbsl
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    std::vector<int> count_vec(nlevs,0);

    SLICEOWNLOOP2
    {
        if(p->flagslice2(i-1,j)<0)
        ++count_vec[p->level];

        if(p->flagslice2(i,j+1)<0)
        ++count_vec[p->level];

        if(p->flagslice2(i,j-1)<0)
        ++count_vec[p->level];

        if(p->flagslice2(i+1,j)<0)
        ++count_vec[p->level];
    }

    p->gcbsl2.resize_levels(nlevs);

    for(int lev=0; lev<nlevs; ++lev)
    p->gcbsl2[lev].assign(count_vec[lev], {});

    // find gcbsl
    count_vec.assign(count_vec.size(), 0);
    SLICEOWNLOOP2
    {
        if(p->flagslice2(i-1,j)<0)
        {
            auto &gcb = p->gcbsl2[p->level][count_vec[p->level]];
            gcb.i=i;
            gcb.j=j;
            gcb.cs=X_NEG;
            GCB_SET_TILE(gcb);
            ++count_vec[p->level];
        }

        if(p->flagslice2(i,j+1)<0)
        {
            auto &gcb = p->gcbsl2[p->level][count_vec[p->level]];
            gcb.i=i;
            gcb.j=j;
            gcb.cs=Y_POS;
            GCB_SET_TILE(gcb);
            ++count_vec[p->level];
        }

        if(p->flagslice2(i,j-1)<0)
        {
            auto &gcb = p->gcbsl2[p->level][count_vec[p->level]];
            gcb.i=i;
            gcb.j=j;
            gcb.cs=Y_NEG;
            GCB_SET_TILE(gcb);
            ++count_vec[p->level];
        }

        if(p->flagslice2(i+1,j)<0)
        {
            auto &gcb = p->gcbsl2[p->level][count_vec[p->level]];
            gcb.i=i;
            gcb.j=j;
            gcb.cs=X_POS;
            GCB_SET_TILE(gcb);
            ++count_vec[p->level];
        }
    }
}
