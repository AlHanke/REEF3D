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
    //int cellcounttopo=0;
    LOOP
    if(a->topo(i,j,k)<p->DZN[KP]) // why is topo for k=1 0.0015 and for k=2 0.03 when topo is located at 0 and k=1 starts at 0 and dz=0.02
    {
        active_topo(i,j,k) = 1.0;
        if(1!=active_box(i,j,k))
        cellcount++;
    }
    //cout<<"Topo cells of part "<<p->mpirank<<": "<<cellcounttopo<<endl;
    
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

    double radius;
    int *tempPosFlag;
    double **tempPos;
    int tempActive=0;
    LOOP
    if(active_topo(i,j,k)>0.0)
    {
        
        for(qn=0;qn<ppcell;++qn)
        if(tempActive<maxparticle)
        {
            radius = p->Q31/2*double(rand() % int(drand/2) + int(drand/2))/drand;
            tempPos[tempActive][0] = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
            tempPos[tempActive][1] = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
            tempPos[tempActive][2] = p->ZN[KP] + p->DZN[KP]*double(rand() % irand)/drand;
            tempPos[tempActive][RADIUS] = radius;
            // verteilen durch ganze zell und alles unter level set weg schmeissen
            // und dann entlang des normalen vectors zur level set ziehen
            // pls_posseed
            // if below flag=0
            if (tempPos[tempActive][2] < p->ZN[KP] + a->topo(i,j,k))
                tempPosFlag[tempActive]=0;
            else
                tempPosFlag[tempActive]=1;
            ++tempActive;
        }
    }
    for (int i=0;i<=tempActive;i++)
        if(1==tempPosFlag[i])
        {
            pos[pactive][0] = tempPos[tempActive][0];
            pos[pactive][1] = tempPos[tempActive][1];
            pos[pactive][2] = tempPos[tempActive][2];
            pos[pactive][RADIUS] = tempPos[tempActive][RADIUS];
            posflag[pactive]=1;
            ++pactive;
        }

    posactive=pactive;
    pcount=pactive;
        // // POS
        //     if(pcount>0)
        //     {
        //         reseeded++;
        //         pcount--;

        //         pos[PC][0] = (double(i) + (rand()%(irand))/drand)*dx;
        //         pos[PC][1] = (double(j) + (rand()%(irand))/drand)*dx;
        //         pos[PC][2] = (double(k) + (rand()%(irand))/drand)*dx;
        //         pos[PC][3] = phipol(p,a,pos[PC][0],pos[PC][1],pos[PC][2]);
        //         posflag[PC]=3;

        //         phival=MAX(((rand()%(irand))/drand)*epsi,rmin);

        //         lambda=1.0;
        //         qq=0;

        //         do
        //         {
        //         normal(a,pos[PC][0],pos[PC][1],pos[PC][2],pos[PC][3]);
        //         pos[PC][0] += lambda*(phival - pos[PC][3])*nvec[0];
        //         pos[PC][1] += lambda*(phival - pos[PC][3])*nvec[1];
        //         pos[PC][2] += lambda*(phival - pos[PC][3])*nvec[2];

        //         ii=int((pos[PC][0])/dx);
        //         jj=int((pos[PC][1])/dx);
        //         kk=int((pos[PC][2])/dx);
        //         check=boundcheck(p,a,ii,jj,kk,0);
        //         if(check==0)
        //         break;

        //         pos[PC][3] = phipol(p,a,pos[PC][0],pos[PC][1],pos[PC][2]);

        //         lambda/=2.0;
        //         ++qq;
        //         }while((pos[PC][3]>epsi || pos[PC][3]<rmin)&& qq<15);
				
		// 		//posradius(p,a,PC);

        //         if((pos[PC][3]>epsi || pos[PC][3]<rmin) || check==0)
        //         {
        //         posflag[PC]=0;
        //         pcount++;
        //         reseeded--;
        //         }
        //     }
			
        //     if(pcount==0 && posactive<maxparticle)
        //     {	
        //         pos[posactive][0] = (double(i)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][1] = (double(j)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][2] = (double(k)  + (rand()%(irand))/drand)*dx;
        //         pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
        //         posflag[posactive]=3;

        //         phival=MAX(((rand()%(irand))/drand)*epsi,rmin);

        //         lambda=1.0;
        //         qq=0;

        //         do
        //         {
        //         normal(a,pos[posactive][0],pos[posactive][1],pos[posactive][2],pos[posactive][3]);
        //         pos[posactive][0] += lambda*(phival - pos[posactive][3])*nvec[0];
        //         pos[posactive][1] += lambda*(phival - pos[posactive][3])*nvec[1];
        //         pos[posactive][2] += lambda*(phival - pos[posactive][3])*nvec[2];

        //         ii=int((pos[posactive][0])/dx);
        //         jj=int((pos[posactive][1])/dx);
        //         kk=int((pos[posactive][2])/dx);
        //         check=boundcheck(p,a,ii,jj,kk,0);
        //         if(check==0)
        //         break;

        //         pos[posactive][3] = phipol(p,a,pos[posactive][0],pos[posactive][1],pos[posactive][2]);
        //         lambda/=2.0;
        //         ++qq;
        //         }while((pos[posactive][3]>epsi || pos[posactive][3]<rmin) && qq<15);


        //         if(pos[posactive][3]<=epsi && pos[posactive][3]>=rmin && check==1)
        //         {
		// 		//posradius(p,a,posactive);
        //         posactive++;
        //         reseeded++;
        //         }
        //     }
			


}


