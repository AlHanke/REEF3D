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

#include"grid_helper.h"
#include"lexer.h"
#include <cassert>

void grid_helper::fillgcb3(lexer *p)
{
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    int q,n;

    ArrayWrapper3D fgc(p);
    fgc.resize(0);

    p->gcb3_count=p->gcb4_count;

    assert(nlevs==1 && "Error: fillgcb3() should only be called when nlevs==1");

    p->gcb3.resize_levels(nlevs);
    for(int lev=0; lev<nlevs; ++lev)
    p->gcb3[lev].assign(p->gcb3_count, {});

    QGCB3
    p->gcb3[p->level][q] = p->gcb4[p->level][q];

    QGC3LOOP
    {
        auto &gcb = p->gcb3[p->level][q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        GCB_APPLY_TILE(gcb, p->level);

        if(gcb.cs==Z_POS)
        fgc(i,j,k)=1;
    }

    QGC3LOOP
    {
        auto &gcb = p->gcb3[p->level][q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        GCB_APPLY_TILE(gcb, p->level);

        if(gcb.cs==Z_POS && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
        gcb.k-=1;
    }

    QGC3LOOP
    {
        auto &gcb = p->gcb3[p->level][q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        GCB_APPLY_TILE(gcb, p->level);

        if(gcb.cs!=Z_POS && fgc(i,j,k)==1 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
        gcb.cs=-abs(gcb.cs);
    }

    GC_TILE_RESET;
}
