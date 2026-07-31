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
#include "reini.h"
#include "6DOF.h"
#include "fdm.h"
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

    #if USE_AMREX
    define_inflow_outflow_ba();
    #endif
}

void lexer::flagini()
{
    control_calc();
	gridsize();

    lexer *p = this;

    flag4.resize(-1);

    p->level = 0;
    TILE_LOOP
    IJKLOOP
    flag4(i,j,k) = flag4_grid[IJK];

    flag4_grid.reset(); // transported into flag4; not needed afterwards

    flag1.resize(OBJ_FLAG);
    flag2.resize(OBJ_FLAG);
    flag3.resize(OBJ_FLAG);
    flag5.resize();

    // boundary conditions
    IO.resize(0);
    IOSL.resize(0);
    DF.resize(1);
#if USE_AMREX
    m_df123 = make_imf(this, 3, &m_df123);
#endif
    // DF1/DF2/DF3: view mode under USE_AMREX (resize() is a no-op; data lives in m_df123)
    DF1.resize();
    DF2.resize();
    DF3.resize();
    Iarray(DFF,imax*jmax*(kmax+1));

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

void lexer::regrid(fdm* a, reini* preini, sixdof* p6dof, ghostcell* pgc, ioflow* pflow)
{
    #if USE_AMREX
    // grid_amrex::regrid_amrex_box_array_and_distribution_mapping(this, a); // Bug with higher levels of static refinement
    // grid_amrex::update_cell_coordinates();
    // grid_amrex::update_cell_spacing();
    // grid_amrex::update_registered_weno(nlevs);
    // grid_amrex::define_inflow_outflow_ba();
    // preini->start(a,this,a->phi,pgc,pflow);
    // lexer* p = this;
    // int counter = 0;
    // PLAINLOOP
    // {
    //     counter++;
    // }
    // veclength += counter - cellnum;
    // cellnum = counter;
    // a->rhsvec.resize(veclength);
    // a->M.resize(veclength);

    #endif
}