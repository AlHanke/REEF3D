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

#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <array>

void lexer::gridini(ghostcell *pgc)
{
    if(G2==1)
    grid::sigma_coord_ini();

    grid::gridspacing(pgc);

    gcd_ini(pgc);

    setup_amrex_geometry(pgc);
}

void lexer::setup_amrex_geometry(ghostcell *pgc)
{
    using namespace amrex;

    Box domain_box(IntVect(AMREX_D_DECL(0, 0, 0)),
                   IntVect(AMREX_D_DECL(gknox-1, gknoy-1, gknoz-1)));

    RealBox physical_domain(AMREX_D_DECL(global_xmin, global_ymin, global_zmin),
                            AMREX_D_DECL(global_xmax, global_ymax, global_zmax));

    std::array<int,AMREX_SPACEDIM> is_periodic = {periodic1, periodic2, periodic3};

    amrex_geometry = Geometry(domain_box, physical_domain, CoordSys::CoordType::cartesian, is_periodic);

    int local_data[6] = {origin_i, origin_j, origin_k, origin_i + knox - 1, origin_j + knoy - 1, origin_k + knoz - 1};

    int all_data[M10][6];
    MPI_Allgather(local_data, 6, MPI_INT, all_data, 6, MPI_INT, pgc->mpi_comm);

    std::vector<Box> all_boxes(M10);
    Vector<int> pmap(M10);
    for (int rank = 0; rank < M10; ++rank) {
        IntVect lo(AMREX_D_DECL(all_data[rank][0], all_data[rank][1], all_data[rank][2]));
        IntVect hi(AMREX_D_DECL(all_data[rank][3], all_data[rank][4], all_data[rank][5]));
        all_boxes[rank] = Box(lo, hi);
        pmap[rank] = rank;
    }

    amrex_box_array = BoxArray(all_boxes.data(), all_boxes.size());
    amrex_distribution_mapping = DistributionMapping(pmap);

    MFIter::allowMultipleMFIters(true);
}

void lexer::flagini()
{
    control_calc();
	gridsize();

	Iarray(flag1,imax*jmax*kmax);
	Iarray(flag2,imax*jmax*kmax);
	Iarray(flag3,imax*jmax*kmax);
    Iarray(flag5,imax*jmax*kmax);
    Iarray(flag,imax*jmax*kmax);

    Iarray(flagsf1,imax*jmax*kmax);
	Iarray(flagsf2,imax*jmax*kmax);
	Iarray(flagsf3,imax*jmax*kmax);
	Iarray(flagsf4,imax*jmax*kmax);

    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
    flagsf1[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf2[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf3[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    flagsf4[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
    }

    // boundary conditions
    Iarray(IO,imax*jmax*kmax);
    Iarray(IOSL,imax*jmax);
    Iarray(DF,imax*jmax*kmax);
    Iarray(DF1,imax*jmax*kmax);
    Iarray(DF2,imax*jmax*kmax);
    Iarray(DF3,imax*jmax*kmax);

    Iarray(DFBED,imax*jmax);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    DFBED[(i-imin)*jmax + j-jmin] = 1;

    // flag
	makeflag(flag1);
	makeflag(flag2);
	makeflag(flag3);

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    IO[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    IOSL[(i-imin)*jmax + j-jmin] = 0;

    for(i=-margin; i<knox+margin; ++i)
    for(j=-margin; j<knoy+margin; ++j)
    for(k=-margin; k<knoz+margin; ++k)
    DF[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin] = 1;

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

void lexer::gcd_ini(ghostcell *pgc)
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
