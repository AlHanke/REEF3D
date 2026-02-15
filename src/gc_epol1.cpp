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

ghostcell::bc_labels ghostcell::gceval1(lexer *p, int gcv, int bc, int cs)
{
    if(gcv==50)
        return bc_labels::NEUMANN;
    //    Velocities

    // Parallel
    //Wall
    else if((bc==21 || bc==22 || (bc==7 && awa_lable==0)) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==10 || gcv==1))
        return gclabel_u;

    else if((bc==21 || bc==22 || (bc==7 && awa_lable==0)) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==110))
        return bc_labels::NOSLIP;

    else if((bc==21 || bc==22 || (bc==7 && awa_lable==0)) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==114))
        return gclabel_u;

    else if((bc==21 || bc==22 || (bc==7 && awa_lable==0)) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==117))
        return bc_labels::NEUMANN;

    // Topo
    else if((bc==5) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==10 || gcv==1))
        return gclabel_utopo;

    else if((bc==5) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==110))
        return bc_labels::NOSLIP;

    else if((bc==5) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==114))
        return gclabel_utopo;

    else if((bc==5) && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==117))
        return gclabel_utopo;

    else if((bc==21 || bc==22 || bc==5) && gcv==14)
        return bc_labels::NEUMANN;

    // Orthogonal
    else if((bc==21 || bc==22 || bc==5) && (cs==1 || cs==4) && (gcv==10 || gcv==1))
        return gclabel_u_orth;

    else if((bc==21 || bc==22 || bc==5) && (cs==1 || cs==4) && gcv==7)
        return bc_labels::NOSLIP;

    //Inflow
    else if((bc==6 && (cs==1 || cs==4) && (gcv==10 || gcv==1 || gcv==7)))
        return gclabel_u_in;

    //Patch
    else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==10 || gcv==1 || gcv==7))
        return bc_labels::NEUMANN;

    //Outflow
    else if((bc==2 && gclabel_outflow==1) && (gcv==10 || gcv==1) && (cs==2 || cs==3 || cs==5 || cs==6))
        return bc_labels::NEUMANN;

    else if((bc==2 && gclabel_outflow==1) && (gcv==10 || gcv==1) && (cs==1 || cs==4))
        return gclabel_u_out;

    else if(((bc==8 || bc==7) && gclabel_outflow==1) && (gcv==10 || gcv==1) && (cs==1 || cs==4) && p->I10==1)
        return bc_labels::NEUMANN;

    //Free Surface
    else if(bc==3 && (cs==2 || cs==3 || cs==5 || cs==6) && (gcv==10 || gcv==17 || gcv==1))
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==1 || cs==4) && (gcv==10 || gcv==17 || gcv==1))
        return gclabel_u_orth;

    else if((bc==9) && cs==6 && (gcv==10 || gcv==17 || gcv==1))
        return bc_labels::NEUMANN;

    // 6DOF
    else if(bc==41 || bc==42 || bc==43)
        return bc_labels::NHPRESS;

    else if(gcv==999)
        return bc_labels::DEBUG;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro1(field& f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    gcdistro(f,ii,jj,kk,dist,gceval1(p,gcv,bc,cs),cs);
}
