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

#include"lexer.h"
#include"ghostcell.h"

#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

void lexer::gridini(ghostcell *pgc)
{
    if(G2==1)
    sigma_coord_ini();

    gridspacing(pgc);
    gcd_ini(pgc);

    setup_amrex_geometry(pgc);
}

void lexer::setup_amrex_geometry(ghostcell *pgc)
{
    using namespace amrex;

    amrex_geometry.resize(nlevs);
    amrex_box_array.resize(nlevs);
    amrex_distribution_mapping.resize(nlevs);
    amrex_bc.resize(nlevs);
    amr_mf.resize(nlevs);

    flag1_imf.resize(nlevs);
    flag2_imf.resize(nlevs);
    flag3_imf.resize(nlevs);
    flag4_imf.resize(nlevs);
    flag7_imf.resize(nlevs);


    int bc_values[6] = {BCType::bogus};
    switch (bcside1)
    {
        case 1:
            bc_values[0] = PhysBCType::inflow;
            break;
        case 3:
            bc_values[0] = BCType::reflect_even;
            break;
        case 6:
            bc_values[0] = BCType::user_1;
            break;
        case 7:
            bc_values[0] = BCType::user_2;
            break;
        case 21: default:
            bc_values[0] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (bcside4)
    {
        case 2:
            bc_values[1] = PhysBCType::outflow;
            break;
        case 3:
            bc_values[1] = BCType::reflect_even;
            break;
        case 6:
            bc_values[1] = BCType::user_1;
            break;
        case 7:
            bc_values[1] = BCType::user_2;
            break;
        case 21: default:
            bc_values[1] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (bcside3)
    {
        case 3:
            bc_values[2] = BCType::reflect_even;
            break;
        case 6:
            bc_values[2] = BCType::user_1;
            break;
        case 7:
            bc_values[2] = BCType::user_2;
            break;
        case 21: default:
            bc_values[2] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (bcside2)
    {
        case 3:
            bc_values[3] = BCType::reflect_even;
            break;
        case 6:
            bc_values[3] = BCType::user_1;
            break;
        case 7:
            bc_values[3] = BCType::user_2;
            break;
        case 21: default:
            bc_values[3] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (bcside5)
    {
        case 3:
            bc_values[4] = BCType::reflect_even;
            break;
        case 21: default:
            bc_values[4] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (bcside6)
    {
        case 3:
            bc_values[5] = BCType::reflect_even;
            break;
        case 21: default:
            bc_values[5] = B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }

    int is_periodic[AMREX_SPACEDIM] = {periodic1, periodic2, periodic3};

    for (int lev = 0; lev < nlevs; lev++)
    {
        RealBox real_box;
        Box domain;
        if (lev == 0)
        {
            real_box.setLo(RealVect(AMREX_D_DECL(global_xmin, global_ymin, global_zmin)));
            real_box.setHi(RealVect(AMREX_D_DECL(global_xmax, global_ymax, global_zmax)));
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(gknox-1, gknoy-1, gknoz-1)));
        }
        else
        {
            // AMR levels not implemented yet
            real_box.setLo(RealVect(AMREX_D_DECL(0, 0, 0)));
            real_box.setHi(RealVect(AMREX_D_DECL(1, 1, 1)));
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(1, 1, 1)));
        }

        amrex_geometry[lev] = Geometry(domain, &real_box, CoordSys::CoordType::cartesian, is_periodic);

        int local_data[6] = {origin_i, origin_j, origin_k, origin_i + knox - 1, origin_j + knoy - 1, origin_k + knoz - 1};

        int all_data[M10][6];
        MPI_Allgather(local_data, 6, MPI_INT, all_data, 6, MPI_INT, pgc->mpi_comm);

        std::vector<Box> all_boxes(M10);
        Vector<int> pmap(M10);
        for (int rank = 0; rank < M10; ++rank)
        {
            IntVect lo(AMREX_D_DECL(all_data[rank][0], all_data[rank][1], all_data[rank][2]));
            IntVect hi(AMREX_D_DECL(all_data[rank][3], all_data[rank][4], all_data[rank][5]));
            all_boxes[rank] = Box(lo, hi);
            pmap[rank] = rank;
        }

        amrex_box_array[lev] = BoxArray(all_boxes.data(), all_boxes.size());
        amrex_distribution_mapping[lev] = DistributionMapping(pmap);

        amrex_bc[lev].resize(n_comp);
        for (int n = 0; n < n_comp; ++n)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (amrex_geometry[lev].Domain().smallEnd(idim) == amrex_geometry[0].Domain().smallEnd(idim))
                {
                    if (amrex_geometry[lev].isPeriodic(idim))
                        amrex_bc[lev][n].setLo(idim, BCType::int_dir);
                    else
                        amrex_bc[lev][n].setLo(idim, bc_values[2*idim]);
                }
                else
                    amrex_bc[lev][n].setLo(idim, BCType::int_dir);

                if (amrex_geometry[lev].Domain().bigEnd(idim) == amrex_geometry[0].Domain().bigEnd(idim))
                {
                    if (amrex_geometry[lev].isPeriodic(idim))
                        amrex_bc[lev][n].setHi(idim, BCType::int_dir);
                    else
                        amrex_bc[lev][n].setHi(idim, bc_values[2*idim + 1]);
                }
                else
                    amrex_bc[lev][n].setHi(idim, BCType::int_dir);

                if (lev == 0 && amrex_geometry[lev].isPeriodic(idim))
                {
                    amrex_bc[lev][n].setLo(idim, BCType::int_dir);
                    amrex_bc[lev][n].setHi(idim, BCType::int_dir);
                }
            }
        }

        amr_mf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 0, margin);

        flag1_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, margin);
        flag2_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, margin);
        flag3_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, margin);
        flag4_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, margin);
        flag7_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, margin);
    }

    amrex::MFIter::allowMultipleMFIters(true);
    level = 0;
    default_mfi = std::make_unique<amrex::MFIter>(amr_mf[level], false);
    amr_mfi = default_mfi.get();
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
