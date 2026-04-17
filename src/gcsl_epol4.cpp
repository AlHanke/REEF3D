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
#include"slice.h"

ghostcell::bc_labels ghostcell::gcsleval4(int gcv, int bc, int cs)
{
    // general Neuman
    if(gcv==1 || gcv==20 || gcv==30 || gcv==40 || gcv==50 || gcv==55)
        return bc_labels::NEUMANN;

    // vertical w
    else if(bc!=INFLOW && (gcv==12 || gcv==24))
        return bc_labels::NEUMANN;

    else if(bc==INFLOW && gcv==12)
        return bc_labels::NOSLIP;

    // pressure 40
    else if((bc==INFLOW || bc==OUTFLOW || bc==SYMMETRY || bc==WAVEGEN || bc==NUMBEACH || bc==WALL || bc==111 || bc==112 || bc==211 || bc==212) && gcv==44)
        return bc_labels::NEUMANN;

    // Potential Ini
    else if((bc==SYMMETRY || bc==WALL) && gcv==49)
        return bc_labels::NEUMANN;

    else if((bc==INFLOW || bc==OUTFLOW || bc==WAVEGEN || bc==NUMBEACH) && gcv==49)
        return bc_labels::POTENTIAL;

    //Patch eta / Hx / Hy
    else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    // eta
    else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54)  && (p->B98<3 || p->A515<=2))
        return bc_labels::NEUMANN;

    else if(bc==OUTFLOW && (gcv==51 || gcv==54))
        return bc_labels::NEUMANN;

    else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else if(bc==NUMBEACH && (gcv==51 || gcv==54 || p->B99==0 || p->B99==1 || ((gcv==52 || gcv==53) && p->B99==3)))
        return bc_labels::NEUMANN;

    else if(bc==NUMBEACH && (gcv==51 || gcv==52 || gcv==53 || gcv==54) && p->B99==4)
        return bc_labels::SOMMERFELD;

    // Fifsf 60 - 3D
    else if(((cs==X_NEG && p->B98<=2) || cs==Y_POS || cs==Y_NEG || (cs==X_POS && p->B99<=2)) && gcv==60)
        return bc_labels::NEUMANN;

    // eta 150
    else if((bc==INFLOW || bc==WAVEGEN) && (gcv==152 || gcv==154))
        return bc_labels::NEUMANN_X;

    else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==151 || gcv==154))
        return bc_labels::NEUMANN_X;

    else if(bc==NUMBEACH && (gcv==151 || gcv==152 || gcv==153 || gcv==154) &&p->B99==3)
        return bc_labels::NEUMANN_X;

    else if(bc==NUMBEACH && (gcv==151 || gcv==152 || gcv==153 || gcv==154) &&p->B99==4)
        return bc_labels::SOMMERFELD;

    else if((bc==SYMMETRY || bc==WALL) && (gcv==151 || gcv==152 || gcv==153 || gcv==154))
        return bc_labels::NEUMANN_X;

    else if(gcv==155)
        return bc_labels::NEUMANN_X;

    // Fifsf 160 - 2D
    else if(((cs==X_NEG && p->B98<=2) || (cs==X_POS && p->B99<=2)) && gcv==160)
        return bc_labels::NEUMANN_X;

    else
        return bc_labels::NONE;
}

void ghostcell::gcsldistro4(lexer *p, slice &f, int ii, int jj, int gcv, int bc, int cs)
{
    gcsldistro(p,f,ii,jj,gcsleval4(gcv,bc,cs),cs);
}


void ghostcell::gcsldistro4int(sliceint &f, int ii, int jj, int cs)
{
    i=ii;
    j=jj;

    gcsl_neumann_int(f,cs);
}

void ghostcell::gcsldistro4Vint(lexer *p, int *f, int ii, int jj, int cs)
{
    i=ii;
    j=jj;

    gcsl_neumann_V_int(p,f,cs);
}
