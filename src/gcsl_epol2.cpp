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

ghostcell::bc_labels ghostcell::gcsleval2(int gcv, int bc, int cs)
{
    // general Neuman
    if(gcv==1 || gcv==40 || gcv==50)
        return bc_labels::NEUMANN;

    //Wall
    // Parallel
    if((bc==INFLOW || bc==WAVEGEN || bc==NUMBEACH || bc==WALL) && (cs==X_NEG || cs==X_POS || cs==Z_NEG || cs==Z_POS) && (gcv==2 || gcv==11 || gcv==21))
        return gclabel_v;

    // Orthogonal
    else if((bc==SYMMETRY || bc==NUMBEACH || bc==WALL) && (cs==Y_POS || cs==Y_NEG) && (gcv==2 || gcv==11 || gcv==21))
        return bc_labels::NOSLIP;

    //Patch
    else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==2 || gcv==8 || gcv==11 || gcv==21))
        return bc_labels::NEUMANN;

    //Outflow
    else if(bc==OUTFLOW && (cs==Y_POS || cs==Y_NEG) && (gcv==2 || gcv==11 || gcv==21))
        return bc_labels::NEUMANN;

    // Symmetry
    else if(bc==SYMMETRY && (cs==X_NEG || cs==X_POS || cs==Z_NEG || cs==Z_POS) && (gcv==2 || gcv==11 || gcv==21))
        return bc_labels::NEUMANN;

    //Hy
    else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54))
        return bc_labels::NEUMANN;

    else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==51 || gcv==54))
        return bc_labels::NEUMANN;

    else if(bc==NUMBEACH && (p->B99==3 || p->B99==4))
        return bc_labels::NEUMANN;

    else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    //Patch Hy
    else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==55 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN_HY;

    else if((bc==222 || bc==212 || bc==122 || bc==112) && (gcv==55 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcsldistro2(lexer *p, slice &f, int ii, int jj, int gcv, int bc, int cs)
{
    gcsldistro(p,f,ii,jj,gcsleval2(gcv,bc,cs),cs);
}

void ghostcell::gcsldistro2int(sliceint &f, int ii, int jj, int cs)
{
    i=ii;
    j=jj;

    gcsl_neumann_int(f,cs);
}
