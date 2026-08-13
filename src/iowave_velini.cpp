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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::velini(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->B98==3 || p->B98==4)
    {
        double area=0.0;
        double Ai=0.0;
        p->Qi=0.0;
        p->Qo=0.0;

        // in
        // LEVEL_LOOP
        GCINLOOP
        {
            i=p->gcin[p->level][n].i;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            GCIN_TILE(n);

            if(a->phi(i-1,j,k)>=0.0)
            {

                if(a->phi(i-1,j,k)>=0.5*p->DXM)
                area=p->DXM*p->DXM;

                if(a->phi(i-1,j,k)<0.5*p->DXM)
                area=p->DXM*(p->DXM*0.5 + a->phi(i-1,j,k));

                Ai+=area;
                p->Qi+=area*a->u(i-1,j,k);
            }
        }
        GC_TILE_RESET;

        double Ui=p->W10/(Ai>1.0e-20?Ai:1.0e20);

        Ui=pgc->globalmax(Ui);


        ULOOP
        a->u(i,j,k)=Ui;

        // LEVEL_LOOP
        GCINLOOP
        {
            i=p->gcin[p->level][n].i;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            GCIN_TILE(n);

            a->u(i-1,j,k)=Ui;
            a->u(i-2,j,k)=Ui;
            a->u(i-3,j,k)=Ui;

            a->v(i-1,j,k)=0.0;
            a->v(i-2,j,k)=0.0;
            a->v(i-3,j,k)=0.0;

            a->w(i-1,j,k)=0.0;
            a->w(i-2,j,k)=0.0;
            a->w(i-3,j,k)=0.0;
        }
        GC_TILE_RESET;
	}
}
