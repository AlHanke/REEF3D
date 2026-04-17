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

ghostcell::bc_labels ghostcell::gcsleval1(int gcv, int bc, int cs)
{
    // general Neuman
    if(gcv==1 || gcv==40 || gcv==50)
        return bc_labels::NEUMANN;
    //Wall
    // Parallel
    if((bc==NUMBEACH || bc==WALL) && (cs==Y_POS || cs==Y_NEG || cs==Z_NEG || cs==Z_POS) && (gcv==1 || gcv==10 || gcv==20))
        return gclabel_u;
    // Orthogonal
    else if((bc==SYMMETRY || (bc==NUMBEACH && !awa_label) || bc==WALL) && (cs==X_NEG || cs==X_POS) && (gcv==1 || gcv==10 || gcv==20))
        return bc_labels::NOSLIP;
    //Outflow
    else if(bc==OUTFLOW && (cs==X_NEG || cs==X_POS) && (gcv==1 || gcv==10 || gcv==20))
        return bc_labels::OUTFLOWBC;
    //Symmetry
    else if(bc==SYMMETRY && (cs==Y_POS || cs==Y_NEG || cs==Z_NEG || cs==Z_POS) && (gcv==1 || gcv==10 || gcv==20))
        return bc_labels::NEUMANN;
    //Hx
    else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54))
        return bc_labels::NEUMANN;

    else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==52 || gcv==53))
        return bc_labels::NEUMANN;

    else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==51 || gcv==54))
        return bc_labels::NEUMANN_HX;

    else if(bc==NUMBEACH && p->B99==3)
        return bc_labels::NEUMANN;

    else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;
    //Patch
    else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==1 || gcv==7 || gcv==10 || gcv==20))
        return bc_labels::NEUMANN;
    //Patch Hx
    else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN_HX;

    else if((bc==222 || bc==212 || bc==122 || bc==112) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcsldistro1(lexer *p, slice &f, int ii, int jj, int gcv, int bc, int cs)
{
    gcsldistro(p,f,ii,jj,gcsleval1(gcv,bc,cs),cs);
}

void ghostcell::gcsldistro1int(sliceint &f, int ii, int jj, int cs)
{
    i=ii;
    j=jj;

    gcsl_neumann_int(f,cs);
}
