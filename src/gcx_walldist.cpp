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
#include"reini.h"
#include"reini_walld.h"

void ghostcell::walldistance(lexer *p, fdm *a, convection *pdisc, reini *preini, ioflow *pflow,  field& walldist)
{
    int ii,jj,kk;
    double xc,yc,zc;
    double xdist,ydist,zdist;

    walldist.setVal(0.0,true);

    walldist.setVal(1.0e9);

    #if USE_AMREX
    walldist.FillBoundary();
    #else
    gcparax(p,walldist,4);
    #endif

    GC4LOOP
    {
        GCB4_TILE(n);

        ii=p->gcb4[p->level][n].i;
        jj=p->gcb4[p->level][n].j;
        kk=p->gcb4[p->level][n].k;

        if( p->gcb4[p->level][n].bc==21)
        {
            if(p->gcb4[p->level][n].cs==1)
            {
                xc = p->XN[IIP];
                yc = p->YP[JJP];
                zc = p->ZP[KKP];

                for(i=ii;i<p->knox;++i)
                {
                    j=jj;
                    k=kk;
                    xdist = fabs(xc - p->XP[IP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),xdist);
                }
            }
            else if(p->gcb4[p->level][n].cs==3)
            {
                xc = p->XP[IIP];
                yc = p->YN[JJP];
                zc = p->ZP[KKP];

                for(j=jj;j<p->knoy;++j)
                {
                    i=ii;
                    k=kk;
                    ydist = fabs(yc - p->YP[JP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),ydist);
                }
            }
            else if(p->gcb4[p->level][n].cs==5)
            {
                xc = p->XP[IIP];
                yc = p->YP[JJP];
                zc = p->ZN[KKP];

                for(k=kk;k<p->knoz;++k)
                {
                    i=ii;
                    j=jj;
                    zdist = fabs(zc - p->ZP[KP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),zdist);
                }
            }
            else if(p->gcb4[p->level][n].cs==4)
            {
                xc = p->XN[IIP1];
                yc = p->YP[JJP];
                zc = p->ZP[KKP];

                for(i=0;i<=ii;++i)
                {
                    j=jj;
                    k=kk;
                    xdist = fabs(xc - p->XP[IP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),xdist);
                }
            }
            else if(p->gcb4[p->level][n].cs==2)
            {
                xc = p->XP[IIP];
                yc = p->YN[JJP1];
                zc = p->ZP[KKP];

                for(j=0;j<=jj;++j)
                {
                    i=ii;
                    k=kk;
                    ydist = fabs(yc - p->YP[JP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),ydist);
                }
            }
            else if(p->gcb4[p->level][n].cs==6)
            {
                xc = p->XP[IIP];
                yc = p->YP[JJP];
                zc = p->ZN[KKP1];

                for(k=0;k<=kk;++k)
                {
                    i=ii;
                    j=jj;
                    zdist = fabs(zc - p->ZP[KP]);
                    PCHECK
                        walldist(i,j,k)=MIN(walldist(i,j,k),zdist);
                }
            }

        }
    }
    GC_TILE_RESET;

    #if USE_AMREX
    walldist.FillBoundary();
    #else
    gcparax(p,walldist,4);
    #endif
    gcparacox(p,walldist);

    // calculate global position of gcb cell
    count=0;
    GC4LOOP
    {
        GCB4_TILE(n);

        if(p->gcb4[p->level][n].bc==21)
        {
            i=p->gcb4[p->level][n].i;
            j=p->gcb4[p->level][n].j;
            k=p->gcb4[p->level][n].k;

            if(p->gcb4[p->level][n].cs==1)
                walldist(i-1,j,k)=-0.5*p->DXN[IP];
            else if(p->gcb4[p->level][n].cs==4)
                walldist(i+1,j,k)=-0.5*p->DXN[IP];
            else if(p->gcb4[p->level][n].cs==3)
                walldist(i,j-1,k)=-0.5*p->DYN[JP];
            else if(p->gcb4[p->level][n].cs==2)
                walldist(i,j+1,k)=-0.5*p->DYN[JP];
            else if(p->gcb4[p->level][n].cs==5)
                walldist(i,j,k-1)=-0.5*p->DZN[KP];
            else if(p->gcb4[p->level][n].cs==6)
                walldist(i,j,k+1)=-0.5*p->DZN[KP];
        }
        else if(p->gcb4[p->level][n].bc==1 || p->gcb4[p->level][n].bc==2 || p->gcb4[p->level][n].bc==3)
        {
            i=p->gcb4[p->level][n].i;
            j=p->gcb4[p->level][n].j;
            k=p->gcb4[p->level][n].k;

            if(p->gcb4[p->level][n].cs==1)
                walldist(i-1,j,k)=walldist(i,j,k);
            else if(p->gcb4[p->level][n].cs==4)
                walldist(i+1,j,k)=walldist(i,j,k);
            else if(p->gcb4[p->level][n].cs==3)
                walldist(i,j-1,k)=walldist(i,j,k);
            else if(p->gcb4[p->level][n].cs==2)
                walldist(i,j+1,k)=walldist(i,j,k);
            else if(p->gcb4[p->level][n].cs==5)
                walldist(i,j,k-1)=walldist(i,j,k);
            else if(p->gcb4[p->level][n].cs==6)
                walldist(i,j,k+1)=walldist(i,j,k);
        }
    }
    GC_TILE_RESET;

    reini_walld reini(p,a);

    reini.start(a,p,walldist,this,pflow);

    #if USE_AMREX
    walldist.FillBoundary();
    #else
    gcparax(p,walldist,4);
    #endif
    gcparacox(p,walldist);
}
