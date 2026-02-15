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

ghostcell::bc_labels ghostcell::gceval2(lexer *p, int gcv, int bc, int cs)
{
    // Velocities
    if(gcv==50)
        return bc_labels::NEUMANN;

    // Parallel
    //Wall
    else if((bc==6 || bc==7 || bc==21) && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && (gcv==11 || gcv==115))
        return gclabel_v;

    else if((bc==6 || bc==7 || bc==21) && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==111)
        return bc_labels::NOSLIP;

    else if(bc==21 && gcv==15)
        return bc_labels::NEUMANN;

    // Orthogonal
    else if((bc==7 || bc==21) && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==11)
        return gclabel_v_orth;

    else if((bc==7 || bc==21) && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==8)
        return bc_labels::NOSLIP;

    //Inflow
    else if(bc==6 && (gcv==8 || gcv==11))
        return gclabel_v_in;

    //Outflow
    else if(bc==2 && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==11 && gclabel_outflow)
        return bc_labels::NEUMANN;

    else if(bc==2 && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==11 && gclabel_outflow)
        return gclabel_v_out;

    //Patch
    else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==8 || gcv==11))
        return bc_labels::NEUMANN;

    //Free Surface
    else if(bc==3 && (cs==dir_labels::X_NEG || cs==dir_labels::X_POS || cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS) && gcv==11)
        return bc_labels::NEUMANN;

    else if(bc==3 && (cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG) && gcv==11)
        return bc_labels::DIRICHLET_ORTH;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro2(field& f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    gcdistro(f,ii,jj,kk,dist,gceval2(p,gcv,bc,cs),cs);
}
