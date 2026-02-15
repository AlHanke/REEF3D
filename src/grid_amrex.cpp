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

    IntVect ref_vec = ref_ratio * IntVect::TheUnitVector();
    if(!j_dir)
        ref_vec[1] = 1;

    for (int lev = 0; lev < nlevs; lev++)
    {
        if (lev == 0) // Overall problem domain
        {
            RealBox real_box;
            real_box.setLo(RealVect(AMREX_D_DECL(p->global_xmin, p->global_ymin, p->global_zmin)));
            real_box.setHi(RealVect(AMREX_D_DECL(p->global_xmax, p->global_ymax, p->global_zmax)));
            Box domain;
            domain.setSmall(IntVect(AMREX_D_DECL(0, 0, 0)));
            domain.setBig(IntVect(AMREX_D_DECL(p->gknox-1, p->gknoy-1, p->gknoz-1)));
            int is_periodic[AMREX_SPACEDIM] = {p->periodic1, p->periodic2, p->periodic3};
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
            if(amrex_box_array[lev].minimalBox() != amrex_geometry[lev].Domain())
            {
                std::cerr << "BoxArray based on grid data does not match the overall problem domain. Check grid parameters." << std::endl;
                exit(1);
            }

            amrex_distribution_mapping[lev] = DistributionMapping(pmap);
        }
        else // Local refinement areas
        {
            std::cerr << "Error: Only single level (no refinement) is currently supported. Exiting." << std::endl;
            exit(1);

            amrex_geometry[lev] = amrex::refine(amrex_geometry[lev-1], ref_vec);

            const double epsion = 1.e-12;
            amrex::BoxList BoxList_lev_n;
            for (const auto &coord_pair : amrex_refined_grid_coords[lev-1])
            {
                // Convert physical coordinates to Level n-1 cell indices
                // Subtract a small epsilon from hi to ensure the index stays within the intended boundary
                amrex::Real lo_phys[3] = {coord_pair.first[0], coord_pair.first[1], coord_pair.first[2]};
                amrex::IntVect lo_idx = amrex_geometry[lev-1].CellIndex(lo_phys);
                amrex::Real hi_phys[3] = {coord_pair.second[0] - epsion, coord_pair.second[1] - epsion, coord_pair.second[2] - epsion};
                amrex::IntVect hi_idx = amrex_geometry[lev-1].CellIndex(hi_phys);
                amrex::Box ref_region_lev_n(lo_idx, hi_idx);
                if(!amrex_box_array[lev-1].contains(ref_region_lev_n))
                {
                    std::cerr << "Refinement region ("<<coord_pair.first[0]<<","<<coord_pair.first[0]<<","<<coord_pair.first[0]<<") to ("<<coord_pair.second[0]<<","<<coord_pair.second[0]<<","<<coord_pair.second[0]<<") at level " << lev << " is not fully contained within the coarser level " << lev-1 << " BoxArray. Check refined grid coordinates." << std::endl;
                    exit(1);
                }
                // Refine the Level n-1 box to Level n index space
                amrex::Box domain_patch = amrex::refine(ref_region_lev_n, ref_vec);
                BoxList_lev_n.push_back(domain_patch);
            }

            amrex_box_array[lev] = amrex::BoxArray(BoxList_lev_n);

            amrex_distribution_mapping[lev] = DistributionMapping(amrex_box_array[lev]);
        }

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

void grid_amrex::update_cell_coordinates()
{
    increment::max_i = (gknox + 2*marge)*pow(ref_ratio,nlevs-1);
    increment::max_j = (gknoy + 2*marge)*pow(ref_ratio,nlevs-1);
    increment::max_k = (gknoz + 2*marge)*pow(ref_ratio,nlevs-1);

    XN.resize((max_i+1)*nlevs,0);
    YN.resize((max_j+1)*nlevs,0);
    ZN.resize((max_k+1)*nlevs,0);

    // Refine the coordinates from level n-1 to level n by inserting midpoints
    for (int lev = 1; lev < nlevs; ++lev)
    {
        for (int i = 0; i < (gknox + 2*marge)*pow(ref_ratio,lev-1) + 1; ++i)
        {
            int idx = 2*i + lev*max_i;
            XN[idx] = XN[i + (lev-1)*max_i];
            idx += 1;
            XN[idx] = XN[i + (lev-1)*max_i]+0.5*(XN[i+1 + (lev-1)*max_i]-XN[i + (lev-1)*max_i]);
        }
        for (int j = 0; j < (gknoy + 2*marge)*pow(ref_ratio,lev-1) +1; ++j)
        {
            int idx = 2*j + lev*max_j;
            YN[idx] = YN[j + (lev-1)*max_j];
            idx += 1;
            YN[idx] = YN[j + (lev-1)*max_j]+0.5*(YN[j+1 + (lev-1)*max_j]-YN[j + (lev-1)*max_j]);
        }
        for (int k = 0; k < (gknoz + 2*marge)*pow(ref_ratio,lev-1) + 1; ++k)
        {
            int idx = 2*k + lev*max_k;
            ZN[idx] = ZN[k + (lev-1)*max_k];
            idx += 1;
            ZN[idx] = ZN[k + (lev-1)*max_k]+0.5*(ZN[k+1 + (lev-1)*max_k]-ZN[k + (lev-1)*max_k]);
        }
    }

    XP.resize((max_i)*nlevs,0);
    YP.resize((max_j)*nlevs,0);
    ZP.resize((max_k)*nlevs,0);

    // Compute cell-centered coordinates as midpoints between node-centered coordinates
    for (int lev = 1; lev < nlevs; ++lev)
    {
        for (int i = 0; i < (gknox + 2*marge)*pow(ref_ratio,lev); ++i)
        {
            int idx = i + lev*max_i;
            XP[idx] = 0.5*(XN[idx]+XN[idx+1]);
        }
        for (int j = 0; j < (gknoy + 2*marge)*pow(ref_ratio,lev); ++j)
        {
            int idx = j + lev*max_j;
            YP[idx] = 0.5*(YN[idx]+YN[idx+1]);
        }
        for (int k = 0; k < (gknoz + 2*marge)*pow(ref_ratio,lev); ++k)
        {
            int idx = k + lev*max_k;
            ZP[idx] = 0.5*(ZN[idx]+ZN[idx+1]);
        }
    }
}

void grid_amrex::update_cell_spacing()
{
    DXN.resize((max_i)*nlevs,0);
    DYN.resize((max_j)*nlevs,0);
    DZN.resize((max_k)*nlevs,0);
    DXP.resize((max_i)*nlevs,0);
    DYP.resize((max_j)*nlevs,0);
    DZP.resize((max_k)*nlevs,0);

    // Compute cell-centered coordinates as midpoints between node-centered coordinates
    for (int lev = 1; lev < nlevs; ++lev)
    {
        for (int i = 0; i < (gknox + 2*marge)*pow(ref_ratio,lev); ++i)
        {
            int idx = i + lev*max_i;
            DXN[idx] = XN[idx+1]-XN[idx];
            DXP[idx] = XP[idx+1]-XP[idx];
        }
        for (int j = 0; j < (gknoy + 2*marge)*pow(ref_ratio,lev); ++j)
        {
            int idx = j + lev*max_j;
            DYN[idx] = YN[idx+1]-YN[idx];
            DYP[idx] = YP[idx+1]-YP[idx];
        }
        for (int k = 0; k < (gknoz + 2*marge)*pow(ref_ratio,lev); ++k)
        {
            int idx = k + lev*max_k;
            DZN[idx] = ZN[idx+1]-ZN[idx];
            DZP[idx] = ZP[idx+1]-ZP[idx];
        }
    }
}
