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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include<math.h>

ghostcell::bc_labels ghostcell::gceval4a(lexer *p, int gcv, int bc, int cs)
{
    // fb
    if(gcv==50)
        return bc_labels::NEUMANN;

    //topo
    else if((bc==3 || bc==5 || bc==6 || bc==7 || bc==8 || bc==21 || bc==22) && (cs==5 || cs==6) && (gcv==151 || gcv==152 || gcv==153))
        return bc_labels::NEUMANN;

    else if((bc==3 || bc==5 || bc==6 || bc==7 || bc==8 || bc==21 || bc==22) && (cs!=5 && cs!=6) && (gcv==151 || gcv==152 || gcv==153))
        return bc_labels::NEUMANN;

    else if(bc==1 && gcv==152)
        return bc_labels::NEUMANN;

    else if(bc==2 && gcv==151)
        return bc_labels::NEUMANN;

    else if(gcv==150 || gcv==154)
     return bc_labels::NEUMANN;

    else if(gcv==159)
        return bc_labels::EXTEND;

    //topo for bedload
    else if((bc==3 || bc==5 || bc==6 || bc==7 || bc==8 || bc==21 || bc==22) && (cs==5 || cs==6) && (gcv==161 || gcv==162 || gcv==163))
     return bc_labels::NEUMANN;

    else if((bc==3 || bc==5 || bc==6 || bc==7 || bc==8 || bc==21 || bc==22) && (cs!=5 && cs!=6) && (gcv==161 || gcv==162 || gcv==163))
     return bc_labels::NEUMANN;

    else if(bc==1 && gcv==162)
        return bc_labels::NEUMANN;

    else if(bc==2 && gcv==161)
        return bc_labels::NEUMANN;

    // Level Set
    else if((bc==5 || bc==6 || bc==7 || bc==8 || bc==9 || bc==21 || bc==22 || bc==41) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else if(bc==3 && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else if(bc==1 && (gcv==52 || gcv==54))
        return bc_labels::NEUMANN;

    else if(bc==2 && (gcv==51 || gcv==54))
        return bc_labels::NEUMANN;

    else if(gcv==50)
        return bc_labels::NEUMANN;

    // porosity
    else if(gcv==1)
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro4a(field& f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    cs = abs(cs);

    gcdistro(f,ii,jj,kk,dist,gceval4a(p,gcv,bc,cs),cs);
}
