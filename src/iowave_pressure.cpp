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
#include"patchBC_interface.h"
#include"heaviside_ls.h"

void iowave::pressure_io(lexer *p, fdm* a, ghostcell *pgc)
{
    pressure_inlet(p,a,pgc);
    pressure_outlet(p,a,pgc);

    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void iowave::pressure_outlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

    // LEVEL_LOOP
    GCOUTLOOP
    {
        i=p->gcout[p->level][n].i;
        j=p->gcout[p->level][n].j;
        k=p->gcout[p->level][n].k;

        GCOUT_TILE(n);

        pval=0.0;

        if(p->B77==1 && p->B99==0)
        {
            if(p->F50==2 || p->F50==3)
            pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
            else if(p->F50==1 || p->F50==4)
            pval=a->press(i,j,k);

            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
        }
        else if(p->B77==1 && (p->B99==1 || p->B99==2))
        {
            pval=a->press(i,j,k);

            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
        }
        else if(p->B77==10)
        {
            double eps,H;

            eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

            H = heaviside_ls(a->phi(i,j,k), eps);

            pval=(1.0-H)*a->press(i,j,k);

            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
        }
    }
    GC_TILE_RESET;
}

void iowave::pressure_inlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

    if(p->B76==0 && p->A10!=5)
    {
        // LEVEL_LOOP
        GCINLOOP
        {
            i=p->gcin[p->level][n].i;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            GCIN_TILE(n);

            if(a->phi(i,j,k)>=0.0)
            pval=(p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
            else
            pval = a->press(i,j,k);

            a->press(i-1,j,k)=pval;
            a->press(i-2,j,k)=pval;
            a->press(i-3,j,k)=pval;
        }
        GC_TILE_RESET;
    }
    else if(p->B76==0 && p->A10==5)
    {
        // LEVEL_LOOP
        GCINLOOP
        {
            i=p->gcin[p->level][n].i;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            GCIN_TILE(n);

            pval = 0.0;

            a->press(i-1,j,k)=pval;
            a->press(i-2,j,k)=pval;
            a->press(i-3,j,k)=pval;
        }
        GC_TILE_RESET;
    }
}

double iowave::local_fsf(lexer *p, fdm *a, ghostcell *pgc)
{
    double wsf=-1.0e20;
    int count=0;

    KLOOP
    PCHECK
    {
        if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
        wsf=MAX(wsf,-(a->phi(i,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());

        if(a->phi(i,j,k)>0.0)
        ++count;
    }

    if(count==0)
    wsf=0.0;

    return wsf;
}
