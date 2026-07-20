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

void grid_helper::fillgcb1(lexer *p)
{
    int q;

    p->Iarray(fgc,imax*jmax*kmax,6);

    p->gcb1_count=p->gcb4_count;
    p->gcb1.assign(p->gcb1_count, {});

    QGCB1
    {
        auto &gcb = p->gcb1[q];

        gcb.i=p->gcb4[q][0];
        gcb.j=p->gcb4[q][1];
        gcb.k=p->gcb4[q][2];
        gcb.cs=p->gcb4[q][3];
        gcb.bc=p->gcb4[q][4];
    }

    QGC1LOOP
    {
        auto &gcb = p->gcb1[q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        fgc[IJK][gcb.cs-1]=1;
    }

    QGC1LOOP
    {
        auto &gcb = p->gcb1[q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        if(gcb.cs==X_POS && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
        gcb.i-=1;
    }

    QGC1LOOP
    {
        auto &gcb = p->gcb1[q];

        i=gcb.i;
        j=gcb.j;
        k=gcb.k;

        if(gcb.cs!=X_POS && fgc[IJK][3]==1 && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
        gcb.cs=-abs(gcb.cs);
    }

    for(int n=0; n<imax*jmax*kmax;n++)
    delete[] fgc[n];

    delete[] fgc;
}
