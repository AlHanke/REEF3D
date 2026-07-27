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

#include"ghostcell.h"
#include"lexer.h"
#include"field.h"

void ghostcell::solid_forcing_lsm(lexer *p, fdm*, field &f)
{
    if(p->X48==1)
    GCDF4LOOP
    {
        GCDF4_TILE(n);
        i=p->gcdf4[p->level][n].i;
        j=p->gcdf4[p->level][n].j;
        k=p->gcdf4[p->level][n].k;

        if(p->gcdf4[p->level][n].cs==1)
        {
            f(i-1,j,k)=f(i,j,k);
            f(i-2,j,k)=f(i,j,k);
        }
        else if(p->gcdf4[p->level][n].cs==4)
        {
            f(i+1,j,k)=f(i,j,k);
            f(i+2,j,k)=f(i,j,k);
        }
        else if(p->gcdf4[p->level][n].cs==3)
        {
            f(i,j-1,k)=f(i,j,k);
            f(i,j-2,k)=f(i,j,k);
        }
        else if(p->gcdf4[p->level][n].cs==2)
        {
            f(i,j+1,k)=f(i,j,k);
            f(i,j+2,k)=f(i,j,k);
        }
        else if(p->gcdf4[p->level][n].cs==5)
        {
            f(i,j,k-1)=f(i,j,k);
            f(i,j,k-2)=f(i,j,k);
        }
        else if(p->gcdf4[p->level][n].cs==6)
        {
            f(i,j,k+1)=f(i,j,k);
            f(i,j,k+2)=f(i,j,k);
        }
    }
    GC_TILE_RESET;
}
