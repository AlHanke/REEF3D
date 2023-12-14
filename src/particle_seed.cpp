/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_f::seed_ini(lexer* p, fdm* a, ghostcell* pgc)
{
    // ini
    LOOP
    {
        active_box(i,j,k) = 0.0;
        active_topo(i,j,k) = 0.0;
    }
    
    // Box
    cellcount=0;
    for(qn=0;qn<p->Q110;++qn)
    LOOP
	if(p->XN[IP]>=p->Q110_xs[qn] && p->XN[IP]<p->Q110_xe[qn]
	&& p->YN[JP]>=p->Q110_ys[qn] && p->YN[JP]<p->Q110_ye[qn]
	&& p->ZN[KP]>=p->Q110_zs[qn] && p->ZN[KP]<p->Q110_ze[qn])
	{
	active_box(i,j,k) = 1.0;
    ++cellcount;
	}

    // Topo
    PLAINLOOP
        if((abs(a->topo(i,j,k))<(p->DZN[KP]*ceil(p->Q102)))&&(a->topo(i,j,k)<=0.25*p->DZN[KP])) //find better comparison to fix numerical drifts
        {
            active_topo(i,j,k) = 1.0;
            if(1!=active_box(i,j,k))
                cellcount++;
        }

    // guess particle demand
    if(p->Q24>0)
        ppcell = p->Q24;
    else
        ppcell = 0;
    
    partnum = cellcount * ppcell;
}

void particle_f::seed(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q110>0)
        posseed_box(p,a,pgc);
	if(p->Q101>0)
        posseed_topo(p,a,pgc);
}


void particle_f::posseed_box(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*(0==p->mpirank?1:p->mpirank));
	
    LOOP
        if(active_box(i,j,k)>0.0)
            for(qn=0;qn<ppcell;++qn)
                if(pactive<maxparticle)
                {
                    pos[pactive][0] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    pos[pactive][1] = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    pos[pactive][2] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
                    pos[pactive][RADIUS] = p->Q31/2*double(rand() % int(drand/2) + int(drand/2))/drand;
                    posflag[pactive]=1;
                    ++pactive;
                }
    
    posactive=pactive;
    pcount=pactive;
}


void particle_f::posseed_topo(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q29>0)
        srand(p->Q29);

    if(p->Q29==0)
        srand((unsigned)time(0)*p->mpirank==0?1:p->mpirank);

    int *tempPosFlag;
    double **tempPos;
    int tempActive=0;
    p->Darray(tempPos,maxparticle,PARTICLE_INFORMATIONS);
    p->Iarray(tempPosFlag,maxparticle);

    PLAINLOOP
        if(active_topo(i,j,k)>0.0)
            for(qn=0;qn<ppcell;++qn)
                if(tempActive<maxparticle)
                {
                    tempPos[tempActive][0] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                    tempPos[tempActive][1] = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                    tempPos[tempActive][2] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
                    tempPos[tempActive][RADIUS] = p->Q31/2*double(rand() % int(drand/2) + int(drand/2))/drand;
                    double ipolTopo = p->ccipol4_b(a->topo,tempPos[tempActive][0],tempPos[tempActive][1],tempPos[tempActive][2]);
                    if (ipolTopo>5e-18||ipolTopo<-p->Q102*p->DZN[KP])
                        tempPosFlag[tempActive]=0;
                    else
                        tempPosFlag[tempActive]=1;
                    ++tempActive;
                }

    for (int i=0;i<tempActive;i++)
        if(1==tempPosFlag[i])
        {
            pos[pactive][0] = tempPos[i][0];
            pos[pactive][1] = tempPos[i][1];
            pos[pactive][2] = tempPos[i][2];
            pos[pactive][RADIUS] = tempPos[i][RADIUS];
            posflag[pactive]=1;
            ++pactive;
        }
    
    posactive=pactive;
    pcount=pactive;

    p->del_Darray(tempPos,maxparticle,PARTICLE_INFORMATIONS);
    p->del_Iarray(tempPosFlag,maxparticle);
}


