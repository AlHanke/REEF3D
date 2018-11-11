/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::rownum4_update(lexer* p, fdm* a,fieldint &rownum4)
{
	p->N4_row=0;
	p->N4_col=0;

    LOOP
	{
    rownum4(i,j,k)=p->N4_row;
    ++p->N4_row;
	++p->N4_col;
	}

    rangex(p,p->range_row4,p->N4_row);

	LOOP
    rownum4(i,j,k)+=p->range_row4[p->mpirank];

}

