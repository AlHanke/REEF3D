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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#if USE_AMREX
#include "fieldint_amrex.h"
#endif

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

    gcd_ini(this,pgc);
    #if USE_AMREX
    define_inflow_outflow_ba();
    #endif
}

void lexer::flagini()
{
    control_calc();
	gridsize();

    flag4.resize();

    lexer* p = this;
    p->level = 0;
    TILE_LOOP
    IJKLOOP
    flag4[IJK] = flag4_grid[IJK];

    flag1.resize(OBJ_FLAG);
    flag2.resize(OBJ_FLAG);
    flag3.resize(OBJ_FLAG);
    flag5.resize();

    flagsf1.resize(1);
    flagsf2.resize(1);
    flagsf3.resize(1);
    flagsf4.resize(1);

    // boundary conditions
    Iarray(IO,imax*jmax*kmax);
    Iarray(IOSL,imax*jmax);
    DF.resize(1);
#if USE_AMREX
    m_df123 = make_imf(this, 3, &m_df123);
#endif
    // DF1/DF2/DF3: view mode under USE_AMREX (resize() is a no-op; data lives in m_df123)
    DF1.resize();
    DF2.resize();
    DF3.resize();

    Iarray(DFBED,imax*jmax);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    DFBED[IJ] = 1;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    IO[IJK] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    IOSL[IJ] = 0;

	x_dir=y_dir=z_dir=1.0;

	if(i_dir==0)
	x_dir=0.0;

	if(j_dir==0)
	y_dir=0.0;

	if(k_dir==0)
	z_dir=0.0;


	if(B98>=3)
	for(n=0;n<gcb4_count;++n)
	if(gcb4[n][4]==6)
	gcb4[n][4]=1;

    // gcdf
    gcdf1_count=gcdf2_count=gcdf3_count=gcdf4_count=1;

    Iarray(gcdf1,gcdf1_count,6);
    Iarray(gcdf2,gcdf2_count,6);
    Iarray(gcdf3,gcdf3_count,6);
    Iarray(gcdf4,gcdf4_count,6);

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

void lexer::gcd_ini(lexer* p, ghostcell *pgc)
{
    for(int q=0;q<gcb4_count;++q)
    {
        if(gcb4[q][3]==1 || gcb4[q][3]==4)
        {
            i=gcb4[q][0];
            gcd4[q] = 0.5*DXP[IP];
        }
        else if(gcb4[q][3]==2 || gcb4[q][3]==3)
        {
            j=gcb4[q][1];
            gcd4[q] = 0.5*DYP[JP];
        }
        else if(gcb4[q][3]==5 || gcb4[q][3]==6)
        {
            k=gcb4[q][2];
            gcd4[q] = 0.5*DZP[KP];
        }
    }
}

void lexer::regrid(fdm* a)
{
    #if USE_AMREX
    grid_amrex::regrid_amrex_box_array_and_distribution_mapping(this, a);
    grid_amrex::update_cell_coordinates();
    grid_amrex::update_cell_spacing();
    grid_amrex::update_registered_weno(nlevs);
    grid_amrex::define_inflow_outflow_ba();
    #endif
}