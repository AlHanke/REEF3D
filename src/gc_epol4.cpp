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
    if(gcv==1 || gcv==30 || gcv==50 || gcv==60 || gcv==150 || gcv==151 || gcv==152 || gcv==153 || gcv==154)
        return bc_labels::NEUMANN;

    // Level Set
    else if((bc==gbc_labels::SYMMETRY || bc==gbc_labels::NUMBEACH || bc==gbc_labels::WALL || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::INFLOW || bc==gbc_labels::WAVEGEN || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==52 || gcv==54))
        return bc_labels::NEUMANN;

    // outflow
    else if((bc==gbc_labels::OUTFLOW || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==51 || (gcv==52 && p->B77==1) || gcv==54))
        return bc_labels::NEUMANN;

    // inflow
    else if((bc==gbc_labels::INFLOW || bc==111 || bc==121 || bc==211 || bc==221) && (gcv==52 || gcv==54))
        return gclabel_lsm_in;

    else if(bc==gbc_labels::WAVEGEN && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return gclabel_lsm_in;

    else if(((bc==gbc_labels::OUTFLOW && !pressout_label) || bc==gbc_labels::SYMMETRY || (bc==gbc_labels::WAVEGEN && !pressin_label) || (bc==gbc_labels::NUMBEACH && !awa_label) ||  bc==gbc_labels::WALL || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return bc_labels::NEUMANN;

    // inflow
    else if(((bc==gbc_labels::INFLOW && !pressin_label) || bc==111 || bc==112 || bc==211 || bc==212) && gcv==40)
        return gclabel_press_in;

    // Turbulence kin
    else if((bc==gbc_labels::OUTFLOW || bc==gbc_labels::WALL) && gcv==20)
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::OUTFLOW || bc==gbc_labels::SYMMETRY) && (cs!=dir_labels::Z_POS || bc!=gbc_labels::SYMMETRY) && gcv==20)
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::SYMMETRY && cs==dir_labels::Z_POS && gcv==20)
        return bc_labels::NOSLIP;

    else if(bc==gbc_labels::SYMMETRY && cs!=dir_labels::Z_POS && gcv==20)
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::INFLOW || bc==gbc_labels::OUTFLOW || bc==gbc_labels::SYMMETRY || bc==gbc_labels::WALL) && gcv==24)
        return bc_labels::NEUMANN;

    else if((bc!=gbc_labels::SYMMETRY || cs!=dir_labels::Z_POS) && gcv==24)
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::WAVEGEN || bc==gbc_labels::NUMBEACH) && gcv==24)
        return bc_labels::NOSLIP;

    // VOF
    else if(bc==gbc_labels::INFLOW && (gcv==72 || gcv==74 || gcv==152))
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::SYMMETRY || bc==gbc_labels::WAVEGEN || bc==gbc_labels::NUMBEACH || bc==gbc_labels::WALL) && (gcv==71 || gcv==72 || gcv==73 || gcv==74 || gcv==151 || gcv==152 || gcv==153))
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::OUTFLOW && (gcv==71 || gcv==74 || gcv==151))
        return bc_labels::NEUMANN;

    // Pk Velocity
    else if(bc==gbc_labels::WALL && (gcv==101 || gcv==102 || gcv==103))
        return bc_labels::NOSLIP;

    // Outflow, Inflow
    else if((bc==gbc_labels::INFLOW || bc==gbc_labels::OUTFLOW || bc==gbc_labels::WAVEGEN || bc==gbc_labels::NUMBEACH) && (gcv==101 || gcv==102 || gcv==103))
        return bc_labels::NEUMANN;

    // Free Surface Uvel
    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==101)
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS) && gcv==101)
        return bc_labels::NOSLIP;

    // Free Surface Vvel
    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==102)
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==102)
        return bc_labels::NOSLIP;

    // Free Surface Wvel
    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::X_NEG || cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG || cs==dir_labels::X_POS) && gcv==103)
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::SYMMETRY && (cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==103)
        return bc_labels::NOSLIP;

    // Heat
    else if(((p->H61==1 && cs==dir_labels::X_NEG) || (p->H62==1 && cs==dir_labels::Y_POS) || (p->H63==1 && cs==dir_labels::Y_NEG) ||
             (p->H64==1 && cs==dir_labels::X_POS) || (p->H65==1 && cs==dir_labels::Z_NEG) || (p->H66==1 && cs==dir_labels::Z_POS)) && gcv==80)
        return bc_labels::HEATBC;

    else if(gcv==80)
        return bc_labels::NEUMANN;

    // Potential Ini
    else if((bc==gbc_labels::WALL || bc==gbc_labels::SYMMETRY) && gcv==49)
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::INFLOW || bc==gbc_labels::OUTFLOW || bc==gbc_labels::WAVEGEN || bc==gbc_labels::NUMBEACH) && gcv==49)
        return bc_labels::POTENTIAL;

    // Potential Waves
    else if((bc==gbc_labels::NUMBEACH || bc==gbc_labels::WALL) && cs!=dir_labels::Z_NEG && gcv==250)
        return bc_labels::NEUMANN;

    else if((bc==gbc_labels::INFLOW || bc==gbc_labels::OUTFLOW || bc==gbc_labels::WAVEGEN || bc==gbc_labels::NUMBEACH) && gcv==250)
        return bc_labels::NEUMANN;

    else if(bc==gbc_labels::SYMMETRY && cs!=dir_labels::Z_POS && gcv==250)
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro4(field &f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    gcdistro(f,ii,jj,kk,dist,gceval4(p,gcv,bc,cs),cs);
}
