/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Alexander Hanke

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "grid_amrex.h"
#include "lexer.h"
#include "ghostcell.h"

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

void grid_amrex::setup_amrex_geometry(lexer* p, ghostcell* pgc)
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
    switch (p->bcside1)
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
            bc_values[0] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (p->bcside4)
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
            bc_values[1] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (p->bcside3)
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
            bc_values[2] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (p->bcside2)
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
            bc_values[3] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (p->bcside5)
    {
        case 3:
            bc_values[4] = BCType::reflect_even;
            break;
        case 21: default:
            bc_values[4] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }
    switch (p->bcside6)
    {
        case 3:
            bc_values[5] = BCType::reflect_even;
            break;
        case 21: default:
            bc_values[5] = p->B20 == 1 ? PhysBCType::slipwall : PhysBCType::noslipwall;
            break;
    }

    int is_periodic[AMREX_SPACEDIM] = {p->periodic1, p->periodic2, p->periodic3};

    for (int lev = 0; lev < nlevs; lev++)
    {
        RealBox real_box;
        Box domain;
        if (lev == 0)
        {
            real_box.setLo(RealVect(AMREX_D_DECL(p->global_xmin, p->global_ymin, p->global_zmin)));
            real_box.setHi(RealVect(AMREX_D_DECL(p->global_xmax, p->global_ymax, p->global_zmax)));
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(p->gknox-1, p->gknoy-1, p->gknoz-1)));
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

        int local_data[6] = {p->origin_i, p->origin_j, p->origin_k, p->origin_i + p->knox - 1, p->origin_j + p->knoy - 1, p->origin_k + p->knoz - 1};

        int all_data[p->M10][6];
        MPI_Allgather(local_data, 6, MPI_INT, all_data, 6, MPI_INT, pgc->mpi_comm);

        std::vector<Box> all_boxes(p->M10);
        Vector<int> pmap(p->M10);
        for (int rank = 0; rank < p->M10; ++rank)
        {
            IntVect lo(AMREX_D_DECL(all_data[rank][0], all_data[rank][1], all_data[rank][2]));
            IntVect hi(AMREX_D_DECL(all_data[rank][3], all_data[rank][4], all_data[rank][5]));
            all_boxes[rank] = Box(lo, hi);
            pmap[rank] = rank;
        }

        amrex_box_array[lev] = BoxArray(all_boxes.data(), all_boxes.size());
        amrex_distribution_mapping[lev] = DistributionMapping(pmap);

        amrex_bc[lev].resize(ncomp);
        for (int n = 0; n < ncomp; ++n)
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

        amr_mf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 0, p->margin);

        flag1_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag2_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag3_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag4_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag7_imf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
    }

    amrex::MFIter::allowMultipleMFIters(true);
    level = 0;
    default_mfi = std::make_unique<amrex::MFIter>(amr_mf[level], false);
    amr_mfi = default_mfi.get();
}