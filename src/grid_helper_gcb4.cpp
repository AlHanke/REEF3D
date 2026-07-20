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
#include"ghostcell.h"
#include"fieldint4.h"

void grid_helper::fillgcb4_wall(lexer *p)
{
    int n;

    fieldint4 cval(p);

    int count=0;

    BASELOOP
    {
        cval(i,j,k)=count;

        ++count;
    }

    GC4LOOP
    {
        GCB4_TILE(n);

        i=p->gcb4[p->level][n].i;
        j=p->gcb4[p->level][n].j;
        k=p->gcb4[p->level][n].k;
        p->gcb4[p->level][n].row=cval(i,j,k);
    }
    GC_TILE_RESET;
}
