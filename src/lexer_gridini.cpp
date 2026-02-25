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
Authors: Hans Bihs, Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#include"lexer.h"
#include"ghostcell.h"

void lexer::gridini(ghostcell *pgc)
{
    #if USE_AMREX
    setup_amrex_geometry(this,pgc);
    #endif

    grid::gridspacing(this, pgc);
    #if USE_AMREX
    if(nlevs > 1)
    {
        if(G2==1)
        {
            std::cerr << "Error: G2=1 is not supported for nlevs > 1" << std::endl;
            pgc->final(true);
        }

        update_cell_coordinates();
        update_cell_spacing();
    }
    #endif

    #if USE_AMREX
    define_inflow_outflow_ba();
    #endif
}

void lexer::flagini()
{
    control_calc();
	gridsize();

    lexer *p = this;

    flag4.resize();

    p->level = 0;
    TILE_LOOP
    IJKLOOP
    flag4[IJK] = flag4_grid[IJK];

    flag4_grid.reset(); // transported into flag4; not needed afterwards

    flag1.resize();
    flag2.resize();
    flag3.resize();
    flag5.resize();

    // boundary conditions
    Iarray(IO,imax*jmax*kmax);
    Iarray(IOSL,imax*jmax);
    DF.resize(1);
    DF1.resize();
    DF2.resize();
    DF3.resize();
    Iarray(DFF,imax*jmax*(kmax+1));

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

    gcsldfeta4.resize(gcsldfeta4_count);
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
