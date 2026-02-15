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

#if USE_AMREX
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
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

void grid_amrex::setup_amrex_geometry(lexer* p, ghostcell* pgc)
{
    bc_type[0] = p->bcside1;
    bc_type[1] = p->bcside4;
    bc_type[2] = p->bcside3;
    bc_type[3] = p->bcside2;
    bc_type[4] = p->bcside5;
    bc_type[5] = p->bcside6;

    using namespace amrex;
    // Initialize AMReX geometry data structures

    amrex_geometry.resize(nlevs);
    amrex_box_array.resize(nlevs);
    amrex_distribution_mapping.resize(nlevs);
    amr_cell_mf.resize(nlevs);

    flag1_iMF.resize(nlevs);
    flag2_iMF.resize(nlevs);
    flag3_iMF.resize(nlevs);
    flag4_iMF.resize(nlevs);
    flag7_iMF.resize(nlevs);

    int is_periodic[AMREX_SPACEDIM] = {p->periodic1, p->periodic2, p->periodic3};

    for (int lev = 0; lev < nlevs; lev++)
    {
        RealBox real_box;
        Box domain;
        if (lev == 0) // Overall problem domain
        {
            real_box.setLo(RealVect(AMREX_D_DECL(p->global_xmin, p->global_ymin, p->global_zmin)));
            real_box.setHi(RealVect(AMREX_D_DECL(p->global_xmax, p->global_ymax, p->global_zmax)));
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(p->gknox-1, p->gknoy-1, p->gknoz-1)));
        }
        else // Local refinement areas
        {
            // AMR levels not implemented yet
            real_box.setLo(RealVect(AMREX_D_DECL(0, 0, 0)));
            real_box.setHi(RealVect(AMREX_D_DECL(1, 1, 1)));
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(1, 1, 1)));
        }

        amrex_geometry[lev] = Geometry(domain, &real_box, CoordSys::CoordType::cartesian, is_periodic);

        // Gather global box information from all ranks to construct BoxArray and DistributionMapping
        int local_data[6] = {p->origin_i, p->origin_j, p->origin_k, p->origin_i + p->knox - 1, p->origin_j + p->knoy - 1, p->origin_k + p->knoz - 1};

        int all_data[p->M10][6];
        MPI_Allgather(local_data, 6, MPI_INT, all_data, 6, MPI_INT, pgc->mpi_comm);

        amrex::Vector<Box> all_boxes(p->M10);
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

        amr_cell_mf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 0, p->margin);

        flag1_iMF[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag2_iMF[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag3_iMF[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag4_iMF[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
        flag7_iMF[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
    }

    amrex::MFIter::allowMultipleMFIters(true);
    level = 0;
    default_cell_mfi = std::make_unique<amrex::MFIter>(amr_cell_mf[level], false);
    amr_cell_mfi = default_cell_mfi.get();

    define_inflow_outflow_ba();
}

void grid_amrex::define_inflow_outflow_ba()
{
    // Intialize inflow and outflow areas
    inflow_ba.resize(nlevs);
    outflow_ba.resize(nlevs);
    inflow_ijk.resize(nlevs);
    outflow_ijk.resize(nlevs);
    for (int lev = 0; lev < nlevs; lev++)
    {
        amrex::BoxList bl_in, bl_out;
        amrex::BoxList original_bl_in, original_bl_out;
        amrex::Vector<int> owner_map_in, owner_map_out;
        const amrex::DistributionMapping& dmap = amrex_distribution_mapping[lev];

        amrex::Box domain = amrex_geometry[lev].Domain();
        for (int i = 0; i < amrex_box_array[lev].size(); ++i)
        {
            amrex::Box b = amrex_box_array[lev][i];
            int owner = dmap[i];

            for (int dir = 0; dir < 3; dir++) // x,y,z directions
            {
                if (b.smallEnd(dir) == domain.smallEnd(dir))
                {
                    amrex::Box b_small = b;
                    b_small.setBig(dir, domain.smallEnd(dir));

                    // Inflow at small end
                    if(bc_type[2*dir]==1 || bc_type[2*dir]==6)
                    {
                        bl_in.push_back(b_small);
                        owner_map_in.push_back(owner);
                        original_bl_in.push_back(b);
                    }

                    // Outflow at small end
                    if (bc_type[2*dir]==2 || bc_type[2*dir]==7)
                    {
                        bl_out.push_back(b_small);
                        owner_map_out.push_back(owner);
                        original_bl_out.push_back(b);
                    }
                }

                if (b.bigEnd(dir) == domain.bigEnd(dir))
                {
                    amrex::Box b_big = b;
                    b_big.setSmall(dir, domain.bigEnd(dir));

                    // Inflow at big end
                    if(bc_type[2*dir+1]==1 || bc_type[2*dir+1]==6)
                    {
                        bl_in.push_back(b_big);
                        owner_map_in.push_back(owner);
                        original_bl_in.push_back(b);
                    }

                    // Outflow at big end
                    if (bc_type[2*dir+1]==2 || bc_type[2*dir+1]==7)
                    {
                        bl_out.push_back(b_big);
                        owner_map_out.push_back(owner);
                        original_bl_out.push_back(b);
                    }
                }
            }
        }

        amrex::BoxArray ba_in(bl_in);
        amrex::DistributionMapping dm_in(owner_map_in);
        inflow_ba[lev].define(ba_in, dm_in, 1, 0);
        inflow_ba[lev].setVal(1);
        inflow_ijk[lev].resize(0);
        for(amrex::MFIter mfi(inflow_ba[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box b = mfi.validbox();
            for (amrex::IntVect iv = b.smallEnd(); iv <= b.bigEnd(); b.next(iv))
            {
                inflow_ijk[lev].push_back(iv-original_bl_in.data()[mfi.index()].smallEnd());
            }
        }
        inflow_ijk[lev].shrink_to_fit();

        amrex::BoxArray ba_out(bl_out);
        amrex::DistributionMapping dm_out(owner_map_out);
        outflow_ba[lev].define(ba_out, dm_out, 1, 0);
        outflow_ba[lev].setVal(1);
        outflow_ijk[lev].resize(0);
        for(amrex::MFIter mfi(outflow_ba[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box b = mfi.validbox();
            for (amrex::IntVect iv = b.smallEnd(); iv <= b.bigEnd(); b.next(iv))
            {
                outflow_ijk[lev].push_back(iv-original_bl_out.data()[mfi.index()].smallEnd());
            }
        }
        outflow_ijk[lev].shrink_to_fit();
    }
}
#endif
