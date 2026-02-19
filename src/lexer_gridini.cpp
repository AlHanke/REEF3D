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

#include"lexer.h"
#include"ghostcell.h"

void lexer::gridini(ghostcell *pgc)
{
    if(G2==1)
    grid::sigma_coord_ini();

    grid::gridspacing(pgc);
}

void lexer::flagini()
{
    control_calc();
	gridsize();

    lexer *p = this;

	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);

    // boundary conditions
    Iarray(IO,imax*jmax*kmax);
    Iarray(IOSL,imax*jmax);
    Iarray(DF,imax*jmax*kmax);
    Iarray(DF1,imax*jmax*kmax);
    Iarray(DF2,imax*jmax*kmax);
    Iarray(DF3,imax*jmax*kmax);

    // flag
	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    IO[IJK] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    IOSL[IJ] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    DF[IJK] = 1;

	if(B98>=3)
	for(n=0;n<gcb4_count;++n)
	if(gcb4[n][4]==6)
	gcb4[n][4]=1;

    // gcdf
    gcdf1_count=gcdf2_count=gcdf3_count=gcdf4_count=1;

    gcdf1.resize(gcdf1_count);
    gcdf2.resize(gcdf2_count);
    gcdf3.resize(gcdf3_count);
    gcdf4.resize(gcdf4_count);

    // gcsldf
    gcsldfeta4_count=1;

    Iarray(gcsldfeta4,gcsldfeta4_count,6);
}

int lexer::conv(double a)
{
	int b,c;
	double d,diff;

	c= int( a);
	d=double(c);
	diff=a-d;

	b=c;

	if(diff>0.5)
	b=c+1;

	if(diff<=-0.5)
	b=c-1;

	return b;
}
