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

void grid_helper::fillgcb1(lexer *p)
{
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    int q;

    ArrayWrapper3D fgc(p);
    fgc.resize(0);

    p->gcb1_count=p->gcb4_count;

    p->gcb1.resize_levels(nlevs);
    LEVEL_LOOP
    {
        LEVEL_DOMAIN_DECL(dom)

        p->gcb1[p->level] = p->gcb4[p->level];

        QGC1LOOP
        {
            auto &gcb = p->gcb1[p->level][q];

            i=gcb.i;
            j=gcb.j;
            k=gcb.k;

            GCB_APPLY_TILE(gcb, p->level);

            if(gcb.cs==X_POS)
            fgc(i,j,k)=1;
        }

        QGC1LOOP
        {
            auto &gcb = p->gcb1[p->level][q];

            i=gcb.i;
            j=gcb.j;
            k=gcb.k;

            GCB_APPLY_TILE(gcb, p->level);

            if(gcb.cs==X_POS && (p->periodic1!=1 || GLOBAL_I<DOMAIN_HI(dom,0)))
            gcb.i-=1;
        }

        QGC1LOOP
        {
            auto &gcb = p->gcb1[p->level][q];

            i=gcb.i;
            j=gcb.j;
            k=gcb.k;

            GCB_APPLY_TILE(gcb, p->level);

            if(gcb.cs!=X_POS && fgc(i,j,k)==1 && (p->periodic1!=1 || GLOBAL_I<DOMAIN_HI(dom,0)))
            gcb.cs=-abs(gcb.cs);
        }
    }
    GC_TILE_RESET;
}