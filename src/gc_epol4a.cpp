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

ghostcell::bc_labels ghostcell::gceval4a(lexer *p, int gcv, int bc)
{
    if(gcv==1 || gcv==50 || gcv==150 || gcv==154)
        return bc_labels::NEUMANN;

    else if((bc==3 || bc==6 || bc==7 || bc==21) && (gcv==151 || gcv==152 || gcv==153))
        return bc_labels::NEUMANN;

    else if(bc==1 && gcv==152)
        return bc_labels::NEUMANN;

    else if(bc==2 && gcv==151)
        return bc_labels::NEUMANN;

    else
        return bc_labels::NONE;
}

void ghostcell::gcdistro4a(field &f, int ii, int jj, int kk, int cs, int bc, double dist, int gcv)
{
    cs = abs(cs);

    gcdistro(f,ii,jj,kk,dist,gceval4a(p,gcv,bc),cs);
}
