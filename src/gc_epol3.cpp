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

ghostcell::bc_labels ghostcell::gceval3(lexer *p, int gcv, int bc, int cs)
{
    // Velocities
    if(gcv==50)
        return bc_labels::NEUMANN;

    // Parallel
    // Wall
    else if((bc==21 || bc==22 || (bc==7 && !awa_label)) && (cs==2 || cs==3 || cs==1 || cs==4) && (gcv==12 || gcv==3))
        return gclabel_w;

    else if((bc==21 || bc==22 || (bc==7 && !awa_label)) && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==112)
        return bc_labels::NOSLIP;

    else if((bc==21 || bc==22 || (bc==7 && !awa_label)) && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==116)
        return gclabel_w;

    else if((bc==21 || bc==22 || (bc==7 && !awa_label)) && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==119)
        return bc_labels::NEUMANN;

    // Topo
    else if(bc==5 && (cs==2 || cs==3 || cs==1 || cs==4) && (gcv==12 || gcv==3))
        return gclabel_wtopo;

    else if(bc==5 && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==112)
        return bc_labels::NOSLIP;

    else if(bc==5 && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==116)
        return gclabel_wtopo;

    else if(bc==5 && (cs==2 || cs==3 || cs==1 || cs==4) && gcv==119)
        return bc_labels::NEUMANN;

    else if((bc==21 || bc==22 || bc==5) && gcv==16)
        return bc_labels::NEUMANN;

    // Othogonal
    else if((bc==21 || bc==22 || bc==5 || (bc==7 && !awa_label)) && (cs==6 || (cs==5 && (p->A10==5 || p->A10==6))) && (gcv==12 || gcv==3))
        return gclabel_w_orth;

    else if((bc==21 || bc==22 || bc==5 || (bc==7 && !awa_label)) && (cs==5 || cs==6) && gcv==9)
        return bc_labels::NOSLIP;

    // Inflow
    else if((bc==6 && (gcv==12 || gcv==3 || gcv==9)))
        return gclabel_w_in;

    // Outflow
    else if((bc==2 && gclabel_outflow==1) && (gcv==12 || gcv==3) && (cs==2 || cs==3 || cs==1 || cs==4))
        return bc_labels::NEUMANN;

    else if((bc==2 && gclabel_outflow==1) && (gcv==12 || gcv==3) && (cs==5 || cs==6))
        return gclabel_w_out;

    // Patch
    else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==12 || gcv==3 || gcv==9))
        return bc_labels::NEUMANN;

    // Free Surface
    else if(bc==3 && (cs==2 || cs==3 || cs==1 || cs==4) && (gcv==12 || gcv==19 || gcv==3))
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==5 || cs==6) && (gcv==12 || gcv==19 || gcv==3) && p->A10!=3 && p->A10!=5)
        return bc_labels::NOSLIP;

    else if(bc==3 && (cs==5 || cs==6) && (gcv==12 || gcv==19 || gcv==3) && p->A10==3)
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==5 || cs==6) && (gcv==12 || gcv==19 || gcv==3 || gcv==112) && p->A10==5)
        return bc_labels::NHPRESS;

    else if(bc==9 && cs==6 && (gcv==12 || gcv==19 || gcv==3))
        return bc_labels::NEUMANN;

    else if(bc==21 && cs==5 && gcv==17)
        return bc_labels::NOSLIP;

    else if((bc!=3 || cs!=6) && gcv==17)
        return bc_labels::NEUMANN;

    // 6DOF
    else if(bc==41 || bc==42 || bc==43)
        return bc_labels::NHPRESS;

    else if(gcv==999)
        return bc_labels::DEBUG;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro3(field &f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    gcdistro(f,ii,jj,kk,dist,gceval3(p,gcv,bc,cs),cs);
}
