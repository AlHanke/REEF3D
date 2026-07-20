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

#include"reini_walld.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reinidisc_f.h"

reini_walld::reini_walld(lexer* p, fdm *a):gradient(p),dab(p)
{
    prdisc = new reinidisc_f(p);
}

reini_walld::~reini_walld()
{
}

void reini_walld::start(fdm* a,lexer* p, field &f, ghostcell* pgc,ioflow* pflow)
{
    pgc->start4(p,f,50);

    int qq;
    QQGC4LOOP
    if(p->gcb4[p->level][qq].bc==21)
    {
        GCB4_TILE(qq);

        i=p->gcb4[p->level][qq].i;
        j=p->gcb4[p->level][qq].j;
        k=p->gcb4[p->level][qq].k;

        if(p->gcb4[p->level][qq].cs==1)
        {
            f(i-1,j,k)=-0.5*p->DXN[IP] - 0.0*p->DXN[IP];
            f(i-2,j,k)=-0.5*p->DXN[IP] - 1.0*p->DXN[IP];
            f(i-3,j,k)=-0.5*p->DXN[IP] - 2.0*p->DXN[IP];
        }
        else if(p->gcb4[p->level][qq].cs==4)
        {
            f(i+1,j,k)=-0.5*p->DXN[IP] - 0.0*p->DXN[IP];
            f(i+2,j,k)=-0.5*p->DXN[IP] - 1.0*p->DXN[IP];
            f(i+3,j,k)=-0.5*p->DXN[IP] - 2.0*p->DXN[IP];
        }
        else if(p->gcb4[p->level][qq].cs==3)
        {
            f(i,j-1,k)=-0.5*p->DYN[JP] - 0.0*p->DYN[JP];
            f(i,j-2,k)=-0.5*p->DYN[JP] - 1.0*p->DYN[JP];
            f(i,j-3,k)=-0.5*p->DYN[JP] - 2.0*p->DYN[JP];
        }
        else if(p->gcb4[p->level][qq].cs==2)
        {
            f(i,j+1,k)=-0.5*p->DYN[JP] - 0.0*p->DYN[JP];
            f(i,j+2,k)=-0.5*p->DYN[JP] - 1.0*p->DYN[JP];
            f(i,j+3,k)=-0.5*p->DYN[JP] - 2.0*p->DYN[JP];
        }
        else if(p->gcb4[p->level][qq].cs==5)
        {
            f(i,j,k-1)=-0.5*p->DZN[KP] - 0.0*p->DZN[KP];
            f(i,j,k-2)=-0.5*p->DZN[KP] - 1.0*p->DZN[KP];
            f(i,j,k-3)=-0.5*p->DZN[KP] - 2.0*p->DZN[KP];
        }
        else if(p->gcb4[p->level][qq].cs==6)
        {
            f(i,j,k+1)=-0.5*p->DZN[KP] - 0.0*p->DZN[KP];
            f(i,j,k+2)=-0.5*p->DZN[KP] - 1.0*p->DZN[KP];
            f(i,j,k+3)=-0.5*p->DZN[KP] - 2.0*p->DZN[KP];
        }
    }
    GC_TILE_RESET;

    QQGC4LOOP
    if(p->gcb4[p->level][qq].bc==1|| p->gcb4[p->level][qq].bc==2|| p->gcb4[p->level][qq].bc==3)
    {
        GCB4_TILE(qq);

        i=p->gcb4[p->level][qq].i;
        j=p->gcb4[p->level][qq].j;
        k=p->gcb4[p->level][qq].k;
        n=p->gcb4[p->level][qq].row;

        if(p->gcb4[p->level][qq].cs==1)
        {
            f(i-1,j,k) = f(i,j,k);
            f(i-2,j,k) = f(i,j,k);
            f(i-3,j,k) = f(i,j,k);
        }
        else if(p->gcb4[p->level][qq].cs==4)
        {
            f(i+1,j,k) = f(i,j,k);
            f(i+2,j,k) = f(i,j,k);
            f(i+3,j,k) = f(i,j,k);
        }
        else if(p->gcb4[p->level][qq].cs==3)
        {
            f(i,j-1,k) = f(i,j,k);
            f(i,j-2,k) = f(i,j,k);
            f(i,j-3,k) = f(i,j,k);
        }
        else if(p->gcb4[p->level][qq].cs==2)
        {
            f(i,j+1,k) = f(i,j,k);
            f(i,j+2,k) = f(i,j,k);
            f(i,j+3,k) = f(i,j,k);
        }
        else if(p->gcb4[p->level][qq].cs==5)
        {
            f(i,j,k-1) = f(i,j,k);
            f(i,j,k-2) = f(i,j,k);
            f(i,j,k-3) = f(i,j,k);
        }
        else if(p->gcb4[p->level][qq].cs==6)
        {
            f(i,j,k+1) = f(i,j,k);
            f(i,j,k+2) = f(i,j,k);
            f(i,j,k+3) = f(i,j,k);
        }
    }
    GC_TILE_RESET;

    double dx;
    dt=1e9;
    if(p->N50==1)
    LOOP
    {
        dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]);

        dt = MIN(dt,0.5*dx);
    }

    reiniter=2*int(p->maxlength/(dt));

    reiniter = pgc->globalimax(reiniter);

    #if USE_AMREX
    f.FillBoundary();
    #else
    pgc->gcparax(p,f,4);
    #endif

    for(int q=0;q<reiniter;++q)
    {

        prdisc->start(p,a,pgc,f,a->L,4);

        if(q==0)
        LOOP
            dab(i,j,k)=a->L(i,j,k);


        LOOP
        {
            f(i,j,k) += dt*0.5*(3.0*a->L(i,j,k) - dab(i,j,k));

            dab(i,j,k)=a->L(i,j,k);
        }

        QQGC4LOOP
        if(p->gcb4[p->level][qq].bc==21)
        {
            GCB4_TILE(qq);

            i=p->gcb4[p->level][qq].i;
            j=p->gcb4[p->level][qq].j;
            k=p->gcb4[p->level][qq].k;
            n=p->gcb4[p->level][qq].row;

            if(p->gcb4[p->level][qq].cs==1 || p->gcb4[p->level][qq].cs==4)
                f(i,j,k) = 0.5*p->DXN[IP];
            else if(p->gcb4[p->level][qq].cs==3 || p->gcb4[p->level][qq].cs==2)
                f(i,j,k) = 0.5*p->DYN[JP];
            else if(p->gcb4[p->level][qq].cs==5 || p->gcb4[p->level][qq].cs==6)
                f(i,j,k) = 0.5*p->DZN[KP];
        }
        GC_TILE_RESET;

        #if USE_AMREX
        f.FillBoundary();
        #else
        pgc->gcparax(p,f,4);
        #endif
    }
}
