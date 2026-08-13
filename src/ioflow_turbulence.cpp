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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void ioflow_f::turbulence_io(lexer *p, fdm* a, ghostcell* pgc)
{
    // LEVEL_LOOP
    GCINLOOP
    {
        i=p->gcin[p->level][n].i;
        j=p->gcin[p->level][n].j;
        k=p->gcin[p->level][n].k;

        GCIN_TILE(n);

        if(a->phi(i-1,j,k)<-1.0*p->F45*p->DXM)
        {
            a->u(i-1,j,k)=a->u(i,j,k);
            a->u(i-2,j,k)=a->u(i,j,k);
            a->u(i-3,j,k)=a->u(i,j,k);

            a->v(i-1,j,k)=a->v(i,j,k);
            a->v(i-2,j,k)=a->v(i,j,k);
            a->v(i-3,j,k)=a->v(i,j,k);

            a->w(i-1,j,k)=a->w(i,j,k);
            a->w(i-2,j,k)=a->w(i,j,k);
            a->w(i-3,j,k)=a->w(i,j,k);
        }
    }
    GC_TILE_RESET;
}
