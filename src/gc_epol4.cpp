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

ghostcell::bc_labels ghostcell::gceval4(lexer *p, int gcv, int bc, int cs)
{
    if(gcv==1 || gcv==50 || gcv==60)
        return bc_labels::NEUMANN;

    // Level Set
    else if((bc==3 || bc==5 || bc==7 || bc==8 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43 || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else if((bc==1 || bc==6 || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==52 || gcv==54))
        return bc_labels::NEUMANN;

    // outflow
    else if((bc==2 || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==51 || (gcv==52 && p->B77==1) || gcv==54))
        return bc_labels::NEUMANN;

    // inflow
    else if((bc==1 || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==52 || gcv==54))
        return gclabel_lsm_in;

    else if(bc==6 && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return gclabel_lsm_in;

    // Pressure
    else if((bc==3 || bc==5 || bc==21 || bc==22 || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return bc_labels::NEUMANN;

    // wavegen
    else if(((bc==6 && !pressin_label) || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return bc_labels::NEUMANN;

    // awa beach
    else if(((bc==7 && !awa_label) || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return bc_labels::NEUMANN;

    // inflow
    else if(((bc==1 && !pressin_label) || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return gclabel_press_in;

    // outflow
    else if(((bc==2 && !pressout_label) || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return bc_labels::NEUMANN;

    // amtosphere
    else if(bc==9 && gcv==40)
        return bc_labels::ATMOSPHERE;

    // Turbulence kin
    else if((bc==5 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && gcv==20)
        return bc_labels::NEUMANN;

    else if((bc==2 || bc==3) && (cs!=dir_labels::Z_POS || bc!=3) && gcv==20)
        return bc_labels::NEUMANN;

    else if(bc==3 && cs==dir_labels::Z_POS && gcv==20)
        return bc_labels::NOSLIP;

    else if((bc==6 || bc==7 || bc==8) && gcv==20)
        return bc_labels::NOSLIP;

    // Turbulence eps
    else if((bc==5 || bc==6 || bc==7 || bc==8 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && gcv==30)
        return bc_labels::NEUMANN;

    else if((bc==2 || bc==3) && gcv==30)
        return bc_labels::NEUMANN;

    else if(bc==1 && gcv==30)
        return bc_labels::NEUMANN;

    // Turbulence eddyv
    else if((bc==5 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && gcv==24)
        return bc_labels::NEUMANN;

    else if((bc==1 || bc==2 || bc==3) && gcv==24)
        return bc_labels::NEUMANN;

    else if(bc==1 && gcv==24)
        return bc_labels::NOSLIP;

    else if(bc==3 && cs==dir_labels::Z_POS && gcv==24)
        return bc_labels::NOSLIP;

    else if((bc!=3 || cs!=dir_labels::Z_POS) && gcv==24)
        return bc_labels::NEUMANN;

    else if((bc==6 || bc==7 || bc==8) && gcv==24)
        return bc_labels::NOSLIP;

    // VOF
    else if((bc==3 || bc==5 || bc==6 || bc==7 || bc==8 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && (gcv==71 || gcv==72 || gcv==73 || gcv==74))
        return bc_labels::NEUMANN;

    else if(bc==1 && (gcv==72 || gcv==74))
        return bc_labels::NEUMANN;

    else if(bc==2 && (gcv==71 || gcv==74))
        return bc_labels::NEUMANN;

    // Pk Velocity
    else if((bc==5 || bc==21 || bc==22 || bc==41) && (gcv==101 || gcv==102 || gcv==103))
        return bc_labels::NOSLIP;

    // Outflow, Inflow
    else if((bc==1 || bc==2 || bc==6 || bc==7 || bc==8) && (gcv==101 || gcv==102 || gcv==103))
        return bc_labels::NEUMANN;

    // Free Surface Uvel
    else if(bc==3 && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==101)
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS) && gcv==101)
        return bc_labels::NOSLIP;

    // Free Surface Vvel
    else if(bc==3 && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==102)
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==102)
        return bc_labels::NOSLIP;

    // Free Surface Wvel
    else if(bc==3 && (cs==dir_labels::X_NEG || cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG || cs==dir_labels::X_POS) && gcv==103)
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==103)
        return bc_labels::NOSLIP;

    // Heat
    else if(((p->H61==1 && cs==dir_labels::X_NEG) || (p->H62==1 && cs==dir_labels::Y_POS) || (p->H63==1 && cs==dir_labels::Y_NEG) ||
             (p->H64==1 && cs==dir_labels::X_POS) || (p->H65==1 && cs==dir_labels::Z_NEG) || (p->H66==1 && cs==dir_labels::Z_POS)) && gcv==80)
        return bc_labels::HEATBC;

    else if(gcv==80)
        return bc_labels::NEUMANN;

    // Potential Ini
    else if((bc==5 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && gcv==49)
        return bc_labels::NEUMANN;

    else if((bc==1 || bc==2 || bc==6 || bc==7 || bc==8) && gcv==49)
        return bc_labels::POTENTIAL;

    else if(bc==3 && gcv==49)
        return bc_labels::NEUMANN;

    // Potential Waves
    else if((bc==5 || bc==7 || bc==8 || bc==9 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43) && cs!=dir_labels::Z_NEG && gcv==250)
        return bc_labels::NEUMANN;

    else if((bc==1 || bc==2 || bc==6 || bc==7) && gcv==250)
        return bc_labels::NEUMANN;

    else if(bc==3 && cs!=dir_labels::Z_POS && gcv==250)
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro4(field &f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    gcdistro(f,ii,jj,kk,dist,gceval4(p,gcv,bc,cs),cs);
}
