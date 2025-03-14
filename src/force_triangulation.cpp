/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more da->solidils.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"force.h"
#include"lexer.h"
#include"fdm.h"

void force::triangulation(lexer *p, fdm* a)
{    
    NDBASELOOP
    {
        eta(i,j,k) = 0.125*(a->solid(i,j,k) + a->solid(i+1,j,k) + a->solid(i,j+1,k) + a->solid(i+1,j+1,k)
                    + a->solid(i,j,k+1) + a->solid(i+1,j,k+1) + a->solid(i,j+1,k+1) + a->solid(i+1,j+1,k+1));
    
        vertice(i,j,k)=-1;

        nodeflag(i,j,k)=0;
    }
    
    BASELOOP
    if(i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke)
    {
        double epsi = interfac*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
        if(fabs(a->solid(i,j,k))<epsi)
        {
            bool check = true;

            if(eta(i,j,k)<0.0 && eta(i-1,j,k)<0.0 && eta(i-1,j-1,k)<0.0 && eta(i,j-1,k)<0.0 &&
               eta(i,j,k-1)<0.0 && eta(i-1,j,k-1)<0.0 && eta(i-1,j-1,k-1)<0.0 && eta(i,j-1,k-1)<0.0)
                check = false;
            
            if(eta(i,j,k)>0.0 && eta(i-1,j,k)>0.0 && eta(i-1,j-1,k)>0.0 && eta(i,j-1,k)>0.0 &&
               eta(i,j,k-1)>0.0 && eta(i-1,j,k-1)>0.0 && eta(i-1,j-1,k-1)>0.0 && eta(i,j-1,k-1)>0.0)
               check = false;

            if(check)
            {
                nodeflag(i,j,k)=1;
                nodeflag(i-1,j,k)=1;
                nodeflag(i-1,j-1,k)=1;
                nodeflag(i,j-1,k)=1;
                nodeflag(i,j,k-1)=1;
                nodeflag(i-1,j,k-1)=1;
                nodeflag(i-1,j-1,k-1)=1;
                nodeflag(i,j-1,k-1)=1;
            }
        }
    }
    
    //--------------------
    int count=0;
    NDBASELOOP
    if(nodeflag(i,j,k)==1)
        ++count;

    numtri = 6*count;
    numvert = count;

    p->Iarray(tri,numtri,4);
    p->Darray(pt,numvert,3);
    p->Darray(ls,numvert);
    p->Iarray(facet,numtri,4);
    p->Iarray(confac,numtri);
    p->Iarray(numfac,numtri);
    p->Iarray(numpt,numtri);
    p->Darray(ccpt,numtri*4,3);

    count=0;
    NDBASELOOP
    if(nodeflag(i,j,k)==1)
    {
        pt[count][0] = p->posnode_x();
        pt[count][1] = p->posnode_y();
        pt[count][2] = p->posnode_z();

        ls[count] = eta(i,j,k);

        vertice(i,j,k) = count;

        ++count;
    }

    // p. 725, 956
    count=0;
    BASELOOP
    if(nodeflag(i,j,k)==1)
    if(nodeflag(i-1,j,k)==1)
    if(nodeflag(i-1,j-1,k)==1)
    if(nodeflag(i,j-1,k)==1)
    if(nodeflag(i,j,k-1)==1)
    if(nodeflag(i-1,j,k-1)==1)
    if(nodeflag(i-1,j-1,k-1)==1)
    if(nodeflag(i,j-1,k-1)==1)
    {
        // 1
        tri[count][0] = vertice(i-1,j-1,k-1);
        tri[count][1] = vertice(i-1,j,k-1);
        tri[count][2] = vertice(i-1,j-1,k);
        tri[count][3] = vertice(i,j-1,k);
        ++count;

        // 2
        tri[count][0] = vertice(i-1,j-1,k-1);
        tri[count][1] = vertice(i,j-1,k-1);
        tri[count][2] = vertice(i-1,j,k-1);
        tri[count][3] = vertice(i,j-1,k);
        ++count;

        // 3
        tri[count][0] = vertice(i-1,j,k-1);
        tri[count][1] = vertice(i,j,k-1);
        tri[count][2] = vertice(i,j-1,k-1);
        tri[count][3] = vertice(i,j-1,k);
        ++count;

        // 4
        tri[count][0] = vertice(i,j,k-1);
        tri[count][1] = vertice(i-1,j,k-1);
        tri[count][2] = vertice(i,j-1,k);
        tri[count][3] = vertice(i,j,k);
        ++count;

        // 5
        tri[count][0] = vertice(i-1,j,k-1);
        tri[count][1] = vertice(i-1,j,k);
        tri[count][2] = vertice(i,j,k);
        tri[count][3] = vertice(i,j-1,k);
        ++count;

        // 6
        tri[count][0] = vertice(i-1,j,k-1);
        tri[count][1] = vertice(i-1,j-1,k);
        tri[count][2] = vertice(i,j-1,k);
        tri[count][3] = vertice(i-1,j,k);
        ++count;
    }
    
    numtri=count;
}
