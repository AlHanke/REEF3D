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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#include "grid_amrex.h"
#include "lexer.h"
#include "ghostcell.h"
#include "fdm.h"
#include "weno3_nug_func.h"
#include "weno_nug_func.h"

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
#include <AMReX_AmrMesh.H>
#include <AMReX_TagBox.H>
#include <AMReX_Array.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

namespace
{
/*!
 * @brief Dynamic AMR mesh that combines static coordinate-based refinement regions
 * (preserving all boxes produced by reef3d_amrmesh_helper) with physics-driven
 * refinement from a level-set field and a Courant-number ceiling derived from the
 * velocity field.
 *
 * Refinement criteria applied in ErrorEst (union — no boxes are ever removed):
 *  1. Static regions: every box specified in @p refined_grid_coords is always tagged.
 *  2. Interface regions: cells where |phi| < @p phi_band_width are tagged.
 *
 * Two independent caps apply to the maximum refinement level:
 *  - @p max_nlevs is a hard ceiling: the effective max level is always
 *    min(max_level, max_nlevs - 1), enforced before any other logic runs.
 *  - Static-region tagging is always applied regardless of either cap.
 */
class reef3d_amrmesh_adaptive final : public amrex::AmrMesh
{
public:
    reef3d_amrmesh_adaptive(
        const amrex::Geometry& level_0_geom,
        const int max_level,
        const amrex::IntVect& level_ref_ratio,
        const amrex::Vector<amrex::Vector<std::pair<amrex::RealVect, amrex::RealVect>>>& refined_grid_coords,
        const amrex::Vector<amrex::MultiFab*>& phi_mf = {},
        const amrex::Vector<amrex::MultiFab*>& fb_mf = {},
        amrex::Real phi_band_width = {},
        amrex::Real fb_band_width = {})
        : amrex::AmrMesh(level_0_geom,
                         make_amr_info(max_level, level_ref_ratio))
        , refined_grid_coords_(refined_grid_coords)
        , phi_mf_(phi_mf)
        , fb_mf_(fb_mf)
        , phi_band_width_(phi_band_width)
        , fb_band_width_(fb_band_width)
    {
    }

private:
    static amrex::AmrInfo make_amr_info(int max_level, const amrex::IntVect& level_ref_ratio)
    {
        amrex::AmrInfo amr_info;
        amr_info.max_level = max_level;
        amr_info.check_input = false;
        amr_info.refine_grid_layout = true;
        amr_info.n_proper = 1;
        amr_info.n_error_buf.assign(max_level + 1, amrex::IntVect::TheZeroVector());
        amr_info.ref_ratio.assign(max_level, level_ref_ratio);
        amr_info.max_grid_size.assign(max_level + 1, amrex::IntVect(AMREX_D_DECL(1048576, 1048576, 1048576)));
        return amr_info;
    }

    void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        tags.setVal(amrex::TagBox::CLEAR);
        if (lev < 0) return;

        // 1. Static refinement regions — identical to reef3d_amrmesh_helper so
        //    that every box it would produce is preserved by this class.
        if (lev < static_cast<int>(refined_grid_coords_.size()))
        {
            constexpr amrex::Real coord_eps = static_cast<amrex::Real>(1.e-12);
            amrex::BoxList refined_regions;

            for (const auto& coord_pair : refined_grid_coords_[lev])
            {
                amrex::Real lo_phys[3] = {coord_pair.first[0], coord_pair.first[1], coord_pair.first[2]};
                amrex::IntVect lo_idx = Geom(lev).CellIndex(lo_phys);

                amrex::Real hi_phys[3] = {coord_pair.second[0] - coord_eps,
                                          coord_pair.second[1] - coord_eps,
                                          coord_pair.second[2] - coord_eps};
                amrex::IntVect hi_idx = Geom(lev).CellIndex(hi_phys);

                amrex::Box tagged_box(lo_idx, hi_idx);
                tagged_box &= Geom(lev).Domain();

                if (tagged_box.ok())
                    refined_regions.push_back(tagged_box);
            }

            if (!refined_regions.isEmpty())
            {
                amrex::BoxArray tag_ba(refined_regions);
                tag_ba.removeOverlap();
                tags.setVal(tag_ba, amrex::TagBox::SET);
            }
        }

        // 2. Fluid-air level-set interface refinement.
        //    Only applied when the resulting finer level (lev+1) is within the
        //    Courant-number ceiling and phi data exists at level 0.
        //    When phi_mf_[lev] has a stale BA (as during MakeNewGrids recursion),
        //    we cascade-interpolate upward from level 0 (whose BA never changes
        //    during regrid) to build phi on the correct BA/DM for this level.
        if (phi_band_width_>0 && !phi_mf_.empty() && phi_mf_[0] != nullptr && !phi_mf_[0]->empty() && lev < static_cast<int>(phi_mf_.size()))
        {
            const amrex::Real band = phi_band_width_;
            const amrex::MultiFab* phi_ptr = nullptr;

            // Direct use when BA matches (first regrid from a stable cycle).
            if (lev < static_cast<int>(phi_mf_.size())
                && phi_mf_[lev] != nullptr
                && !phi_mf_[lev]->empty()
                && phi_mf_[lev]->boxArray() == boxArray(lev))
            {
                phi_ptr = phi_mf_[lev];
            }

            // Cascade interpolation from level 0 when phi at this level is
            // unavailable or on a stale BA. reserve() prevents reallocation so
            // cur_coarse pointer stays valid across iterations.
            std::vector<amrex::MultiFab> interp_chain;
            if (phi_ptr == nullptr)
            {
                amrex::Vector<amrex::BCRec> bcrecs(1);
                for (auto& bc : bcrecs)
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    { bc.setLo(d, amrex::BCType::foextrap); bc.setHi(d, amrex::BCType::foextrap); }
                amrex::PhysBCFunctNoOp null_bc;

                interp_chain.reserve(lev);
                const amrex::MultiFab* cur_coarse = phi_mf_[0];

                for (int cl = 0; cl < lev; ++cl)
                {
                    // Final step lands on boxArray(lev)/DistributionMap(lev) so
                    // the resulting MF has the same BA/DM as tags — no mismatch.
                    amrex::BoxArray fine_ba = (cl + 1 == lev)
                        ? boxArray(lev)
                        : amrex::refine(cur_coarse->boxArray(), refRatio(cl));
                    amrex::DistributionMapping fine_dm = (cl + 1 == lev)
                        ? DistributionMap(lev)
                        : amrex::DistributionMapping(fine_ba);

                    interp_chain.emplace_back(fine_ba, fine_dm, 1, 0);
                    amrex::InterpFromCoarseLevel(
                        interp_chain.back(), 0.0, *cur_coarse, 0, 0, 1,
                        Geom(cl), Geom(cl + 1),
                        null_bc, 0, null_bc, 0,
                        refRatio(cl), &amrex::cell_cons_interp, bcrecs, 0);

                    // Copy existing phi data at this intermediate level on top of
                    // the interpolated values so that known fine-level data is used.
                    const int fine_lev = cl + 1;
                    if (fine_lev < static_cast<int>(phi_mf_.size())
                        && phi_mf_[fine_lev] != nullptr
                        && !phi_mf_[fine_lev]->empty())
                    {
                        interp_chain.back().ParallelCopy(
                            *phi_mf_[fine_lev], 0, 0, 1,
                            phi_mf_[fine_lev]->nGrowVect(),
                            amrex::IntVect::TheZeroVector());
                    }

                    cur_coarse = &interp_chain.back();
                }
                phi_ptr = cur_coarse;
            }

            if (phi_ptr != nullptr)
            {
                for (amrex::MFIter mfi(*phi_ptr); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    const auto phi_arr   = phi_ptr->const_array(mfi);
                    auto tag_arr         = tags.array(mfi);

                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            if (std::abs(phi_arr(i, j, k, 0)) < band)
                                tag_arr(i, j, k) = amrex::TagBox::SET;
                        });
                }
            }
        }

        // 3. Fluid-6DOF level-set interface refinement.
        //    Only applied when the resulting finer level (lev+1) is within the
        //    Courant-number ceiling and phi data exists at level 0.
        //    When phi_mf_[lev] has a stale BA (as during MakeNewGrids recursion),
        //    we cascade-interpolate upward from level 0 (whose BA never changes
        //    during regrid) to build phi on the correct BA/DM for this level.
        if (fb_band_width_>0 && !fb_mf_.empty() && fb_mf_[0] != nullptr && !fb_mf_[0]->empty() && lev < static_cast<int>(fb_mf_.size()) && lev == 0)
        {
            const amrex::Real band = fb_band_width_;
            const amrex::MultiFab* fb_ptr = nullptr;

            // Direct use when BA matches (first regrid from a stable cycle).
            if (lev < static_cast<int>(fb_mf_.size())
                && fb_mf_[lev] != nullptr
                && !fb_mf_[lev]->empty()
                && fb_mf_[lev]->boxArray() == boxArray(lev))
            {
                fb_ptr = fb_mf_[lev];
            }

            // Cascade interpolation from level 0 when fb at this level is
            // unavailable or on a stale BA. reserve() prevents reallocation so
            // cur_coarse pointer stays valid across iterations.
            std::vector<amrex::MultiFab> interp_chain;
            if (fb_ptr == nullptr)
            {
                amrex::Vector<amrex::BCRec> bcrecs(1);
                for (auto& bc : bcrecs)
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    { bc.setLo(d, amrex::BCType::foextrap); bc.setHi(d, amrex::BCType::foextrap); }
                amrex::PhysBCFunctNoOp null_bc;

                interp_chain.reserve(lev);
                const amrex::MultiFab* cur_coarse = fb_mf_[0];

                for (int cl = 0; cl < lev; ++cl)
                {
                    // Final step lands on boxArray(lev)/DistributionMap(lev) so
                    // the resulting MF has the same BA/DM as tags — no mismatch.
                    amrex::BoxArray fine_ba = (cl + 1 == lev)
                        ? boxArray(lev)
                        : amrex::refine(cur_coarse->boxArray(), refRatio(cl));
                    amrex::DistributionMapping fine_dm = (cl + 1 == lev)
                        ? DistributionMap(lev)
                        : amrex::DistributionMapping(fine_ba);

                    interp_chain.emplace_back(fine_ba, fine_dm, 1, 0);
                    amrex::InterpFromCoarseLevel(
                        interp_chain.back(), 0.0, *cur_coarse, 0, 0, 1,
                        Geom(cl), Geom(cl + 1),
                        null_bc, 0, null_bc, 0,
                        refRatio(cl), &amrex::cell_cons_interp, bcrecs, 0);

                    // Copy existing fb data at this intermediate level on top of
                    // the interpolated values so that known fine-level data is used.
                    const int fine_lev = cl + 1;
                    if (fine_lev < static_cast<int>(fb_mf_.size())
                        && fb_mf_[fine_lev] != nullptr
                        && !fb_mf_[fine_lev]->empty())
                    {
                        interp_chain.back().ParallelCopy(
                            *fb_mf_[fine_lev], 0, 0, 1,
                            fb_mf_[fine_lev]->nGrowVect(),
                            amrex::IntVect::TheZeroVector());
                    }

                    cur_coarse = &interp_chain.back();
                }
                fb_ptr = cur_coarse;
            }

            if (fb_ptr != nullptr)
            {
                for (amrex::MFIter mfi(*fb_ptr); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    const auto fb_arr   = fb_ptr->const_array(mfi);
                    auto tag_arr         = tags.array(mfi);

                    amrex::ParallelFor(bx,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            if (std::abs(fb_arr(i, j, k, 0)) < band)
                                tag_arr(i, j, k) = amrex::TagBox::SET;
                        });
                }
            }
        }
    }

    const amrex::Vector<amrex::Vector<std::pair<amrex::RealVect, amrex::RealVect>>>& refined_grid_coords_;
    const amrex::Vector<amrex::MultiFab*>& phi_mf_;
    const amrex::Vector<amrex::MultiFab*>& fb_mf_;
    const amrex::Real phi_band_width_;
    const amrex::Real fb_band_width_;
};
}

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

    ref_vec = ref_ratio * IntVect::TheUnitVector();
    if(!i_dir)
        ref_vec[0] = 1;
    if(!j_dir)
        ref_vec[1] = 1;
    if(!k_dir)
        ref_vec[2] = 1;

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
        else
        {
            amrex_geometry[lev] = amrex::refine(amrex_geometry[lev-1], ref_vec);
        }
    }

    create_amrex_box_array_and_distribution_mapping_level_n();

    for (int lev = 0; lev < nlevs; lev++)
    {
        amr_cell_mf[lev].define(amrex_box_array[lev], amrex_distribution_mapping[lev], 1, p->margin);
    }

    const auto ghost_vec = amrex::IntVect(p->margin, p->margin, p->margin);
    for (int lev = nlevs-1; lev >= 0; lev--)
    {
        if(lev == nlevs-1)
            amr_cell_mf[lev].setVal(0);
        else
            amr_cell_mf[lev] = amrex::makeFineMask(amr_cell_mf[lev], amr_cell_mf[lev+1], ghost_vec, ref_vec, amrex_geometry[lev].periodicity(), 1, 0);
    }

    output_amrex_level_info();

    amrex::MFIter::allowMultipleMFIters(true);
    level = 0;
    default_cell_mfi = std::make_unique<amrex::MFIter>(amr_cell_mf[level], false);
    set_tile_mfi(default_cell_mfi.get());
}

void grid_amrex::output_amrex_level_info()
{
    if (amrex::ParallelDescriptor::MyProc() == 0)
    {
        std::cout<<"Number of AMR levels: "<<nlevs<<std::endl;
    }

    for (int lev = 0; lev < nlevs; ++lev)
    {
        if (amrex::ParallelDescriptor::MyProc() == 0)
        {
            std::cout << "AMReX level " << lev << " cells: " << amrex_box_array[lev].numPts() << std::endl;
            // #ifndef NDEBUG
            // std::cout << "BoxArray: " << amrex_box_array[lev] << std::endl;
            // std::cout << "DistributionMapping: " << amrex_distribution_mapping[lev] << std::endl;
            // #endif
        }
    }
}

/*!
 * @brief Update the BoxArrays and DistributionMappings for finer leves
 * Ensures that all finer level BoxArrays are properly nested within the coarser level BoxArrays.
*/
void grid_amrex::create_amrex_box_array_and_distribution_mapping_level_n()
{
    if (nlevs <= 1)
        return;

    if (static_cast<int>(amrex_refined_grid_coords.size()) < nlevs - 1)
    {
        std::cerr << "Insufficient amrex_refined_grid_coords entries for " << nlevs << " AMR levels." << std::endl;
        exit(1);
    }

    reef3d_amrmesh_adaptive amr_mesh_helper(amrex_geometry[0], nlevs - 1, ref_vec, amrex_refined_grid_coords);
    amr_mesh_helper.SetGeometry(0, amrex_geometry[0]);
    amr_mesh_helper.SetBoxArray(0, amrex_box_array[0]);
    amr_mesh_helper.SetDistributionMap(0, amrex_distribution_mapping[0]);
    amr_mesh_helper.SetFinestLevel(0);

    for (int lev = 1; lev < nlevs; ++lev)
    {
        amr_mesh_helper.SetGeometry(lev, amrex_geometry[lev]);

        amrex::Vector<amrex::BoxArray> new_grids(nlevs);
        int new_finest = lev - 1;
        amr_mesh_helper.MakeNewGrids(lev - 1, 0.0, new_finest, new_grids);

        if (new_finest < lev || new_grids[lev].empty())
        {
            std::cerr << "Failed to build AMR level " << lev << " from amrex_refined_grid_coords using AMReX proper nesting." << std::endl;
            exit(1);
        }

        amrex_box_array[lev] = new_grids[lev];

        const int nprocs = amrex::ParallelDescriptor::NProcs();
        if (nprocs > 1)
            amr_mesh_helper.ChopGrids(lev, amrex_box_array[lev], nprocs);

        amrex_distribution_mapping[lev] = amr_mesh_helper.MakeDistributionMap(lev, amrex_box_array[lev]);

        amr_mesh_helper.SetBoxArray(lev, amrex_box_array[lev]);
        amr_mesh_helper.SetDistributionMap(lev, amrex_distribution_mapping[lev]);
        amr_mesh_helper.SetFinestLevel(lev);
    }
}

/*!
 * @brief Adaptive regrid: update BoxArrays and DistributionMappings for all finer
 * levels using physics-driven refinement.
 *
 * Mirrors create_amrex_box_array_and_distribution_mapping_level_n() but drives
 * ErrorEst through reef3d_amrmesh_adaptive, which adds interface-based tagging
 * (|phi| < p->psi) on top of the static coordinate regions, subject to a
 * Courant-number ceiling computed from the current velocity field and p->dt.
 *
 * Level 0 BoxArray / DistributionMapping and all Geometry objects are assumed
 * to be already set up (e.g. by setup_amrex_geometry).
 *
 * @param p  lexer — provides dt, psi (interface half-width), max_nlevs, ref_vec.
 * @param a  fdm   — provides phi, u, v, w MultiFabs at each AMR level.
 */
void grid_amrex::regrid_amrex_box_array_and_distribution_mapping(lexer* p, fdm* a)
{
    const int old_nlevs = nlevs;

    // Build non-owning pointer vectors from CURRENT (old) levels only.
    // For any new level that gets added, phi_mf_[lev] will be null — the
    // null-check in reef3d_amrmesh_adaptive::ErrorEst skips phi tagging safely.
    amrex::Vector<amrex::MultiFab*> phi_mfs(old_nlevs);
    amrex::Vector<amrex::MultiFab*> fb_mfs(old_nlevs);
    for (int lev = 0; lev < old_nlevs; ++lev)
    {
        phi_mfs[lev] = &a->phi.GetMultiFab(lev);
        fb_mfs[lev] = &a->fb.GetMultiFab(lev);
    }

    i=j=k=0;
    // Construct with max_nlevs-1 (not nlevs-1) so the mesh can detect that
    // additional levels are needed, not just reduce the existing count.
    reef3d_amrmesh_adaptive amr_mesh_adaptive(
        amrex_geometry[0],
        max_nlevs - 1,
        ref_vec,
        amrex_refined_grid_coords,
        phi_mfs,
        fb_mfs,
        static_cast<amrex::Real>(0/*p->F45*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP])*/),
        static_cast<amrex::Real>(0/*1.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP])*/));

    // Extend geometry/BA/DM vectors to full max_nlevs so MakeNewGrids can
    // recurse through all levels without hitting unregistered geometry.
    // Limit probing to avoid container-overflow issues in fine-level error estimation.
    const int max_probe = std::min(2, max_nlevs);  // Limit to level 1 (avoid level 2+ probing)
    if(p->mpirank == 0) std::cout<< "Probing up to " << max_probe << std::endl;

    amrex_geometry.resize(max_nlevs);
    for (int lev = old_nlevs; lev < max_nlevs; ++lev)
        amrex_geometry[lev] = amrex::refine(amrex_geometry[lev-1], ref_vec);

    // Initialize BoxArray and DistributionMapping vectors properly.
    // Ensure level 0 entries exist with valid (possibly empty) values.
    if (amrex_box_array.empty())
        amrex_box_array.push_back(amrex::BoxArray());
    if (amrex_distribution_mapping.empty())
        amrex_distribution_mapping.push_back(amrex::DistributionMapping());

    // Register level 0 (BoxArray / DistributionMapping are unchanged).
    amr_mesh_adaptive.SetGeometry(0, amrex_geometry[0]);
    amr_mesh_adaptive.SetBoxArray(0, amrex_box_array[0]);
    amr_mesh_adaptive.SetDistributionMap(0, amrex_distribution_mapping[0]);
    amr_mesh_adaptive.SetFinestLevel(0);

    // Create a dummy BoxArray and DistributionMapping for probe-only levels.
    // This allows ErrorEst to safely reference higher levels during MakeNewGrids
    // without triggering AMReX container-overflow errors.
    amrex::BoxArray dummy_ba(amrex::Box::TheUnitBox());
    amrex::DistributionMapping dummy_dm(dummy_ba);

    for (int lev = 1; lev < max_nlevs; ++lev)
    {
        if ((int)amrex_box_array.size() <= lev)
            amrex_box_array.push_back(dummy_ba);
        if ((int)amrex_distribution_mapping.size() <= lev)
            amrex_distribution_mapping.push_back(dummy_dm);
    }

    // Pre-allocate all registered MultiFabs to max_nlevs with dummy BoxArrays
    // so ErrorEst callbacks can safely access them during probing.
    for (auto& e : mf_registry)
    {
        while ((int)e.mf->size() < max_nlevs)
            e.mf->emplace_back(dummy_ba, dummy_dm, e.ncomp, 0);
    }
    for (auto& e : imf_registry)
    {
        while ((int)e.mf->size() < max_nlevs)
            e.mf->emplace_back(dummy_ba, dummy_dm, e.ncomp, 0);
    }

    // Register all geometries upfront so MakeNewGrids can recurse to max_level
    // safely (it calls Geom(lev) for domain-clipping at each recursion depth).
    for (int lev = 1; lev < max_nlevs; ++lev)
        amr_mesh_adaptive.SetGeometry(lev, amrex_geometry[lev]);

    // Probe each candidate level; stop as soon as no cells are tagged.
    int actual_finest = 0;
    const int nprocs  = amrex::ParallelDescriptor::NProcs();
    for (int lev = 1; lev < max_probe; ++lev)
    {
        if(p->mpirank == 0) std::cout << "Level " << lev << std::endl;

        // Ensure vectors have space for this level before accessing.
        while ((int)amrex_box_array.size() <= lev)
            amrex_box_array.push_back(amrex::BoxArray());
        while ((int)amrex_distribution_mapping.size() <= lev)
            amrex_distribution_mapping.push_back(amrex::DistributionMapping());

        // Skip probing if the refined grid domain would exceed a practical limit.
        // This prevents overflow issues in AMReX's error estimation at very fine levels.
        amrex::IntVect domain_cells = amrex_geometry[lev].Domain().length();
        if (domain_cells[0] <= 2 || domain_cells[1] <= 2 || domain_cells[2] <= 2)
            break;  // Grid is too fine to subdivide further

        amrex::Vector<amrex::BoxArray> new_grids(max_nlevs);
        int new_finest = lev - 1;
        amr_mesh_adaptive.MakeNewGrids(lev - 1, static_cast<amrex::Real>(p->simtime),
                                       new_finest, new_grids);

        if (new_finest < lev || new_grids[lev].empty())
            break;  // no tagging at this level — no finer levels either

        // Compare the NEW (post-chop) BoxArray against the existing one. Save the old
        // BoxArray first so the comparison is against what we actually had, not the
        // value we are about to overwrite.
        const amrex::BoxArray old_ba = amrex_box_array[lev];

        amrex_box_array[lev] = new_grids[lev];

        if (nprocs > 1)
            amr_mesh_adaptive.ChopGrids(lev, amrex_box_array[lev], nprocs);

        const bool ba_same = (old_ba == amrex_box_array[lev]);
        if(!ba_same) changed = true;

        // DIAGNOSTIC (REEF_FORCE_RANK_FLIP): force a PURE rank change with the BoxArray
        // unchanged, to isolate whether a rank flip alone breaks the C-F projection.
        // Cyclically shift every box's owner by +1 (mod nprocs) and install that DM even
        // when ba_same, so the patch lands on a different rank with identical geometry.
        // Run with REEF_PROJ_CHECK=1: if projcheck jumps to ~4e-3 the C-F path has a
        // genuine cross-rank bug (ordering/staleness); if it stays ~1e-19 the rank change
        // is transparent and the dm-preservation fix merely avoids rebuild churn.
        if(ba_same && std::getenv("REEF_FORCE_RANK_FLIP") && nprocs > 1)
        {
            amrex::Vector<int> pmap(amrex_box_array[lev].size());
            for(int b = 0; b < (int)pmap.size(); ++b)
                pmap[b] = (amrex_distribution_mapping[lev][b] + 1) % nprocs;
            amrex_distribution_mapping[lev] = amrex::DistributionMapping(pmap);
            if(p->mpirank == 0)
                std::cout << "  [rankflip] lev " << lev << " forced +1 owner shift" << std::endl;
            changed = true;
        }
        // Preserve the existing DistributionMapping when the BoxArray is unchanged. The
        // fine patch is a single box owned by one rank; MakeDistributionMap reassigns that
        // box (knapsack/SFC over box index) and can move it to a DIFFERENT rank even when
        // the box itself is identical. That flip changes which rank owns the fine side of
        // every C-F link, and the composite-solve C-F coupling then assembles inconsistently
        // across the new coarse/fine rank pairing -> projcheck residual ~4e-3 at the C-F ->
        // the projection injects velocity from rest. Only rebuild the DM when the boxes
        // actually changed (or the current DM no longer matches the BoxArray size).
        else if(!ba_same
           || amrex_distribution_mapping[lev].size() != amrex_box_array[lev].size())
        {
            amrex_distribution_mapping[lev] =
                amr_mesh_adaptive.MakeDistributionMap(lev, amrex_box_array[lev]);
            changed = true;
        }

        amr_mesh_adaptive.SetBoxArray(lev, amrex_box_array[lev]);
        amr_mesh_adaptive.SetDistributionMap(lev, amrex_distribution_mapping[lev]);
        amr_mesh_adaptive.SetFinestLevel(lev);
        actual_finest = lev;
        if(p->mpirank == 0) std::cout << "Level " << lev << " tagged cells: " << amrex_box_array[lev].numPts() << std::endl;
    }

    const int new_nlevs = actual_finest + 1;

    if (new_nlevs != old_nlevs && p->mpirank == 0)
        std::cout << "AMReX adaptive regrid: levels " << old_nlevs
                  << " -> " << new_nlevs << std::endl;

    // Clear any unused levels above new_nlevs to avoid accessing invalid BoxArrays/DMs.
    // This must be done BEFORE trimming vectors so AMReX cleanup can happen properly.
    for (int lev = new_nlevs; lev < (int)amrex_box_array.size(); ++lev)
    {
        amrex_box_array[lev] = amrex::BoxArray();
        amrex_distribution_mapping[lev] = amrex::DistributionMapping();
    }

    // Trim vectors to only include levels that were actually created by tagging.
    amrex_geometry.resize(new_nlevs);
    amrex_box_array.resize(new_nlevs);
    amrex_distribution_mapping.resize(new_nlevs);

    // Trim registered MultiFabs back to actual level count to avoid dangling dummy MFs.
    for (auto& e : mf_registry)
        e.mf->resize(new_nlevs);
    for (auto& e : imf_registry)
        e.mf->resize(new_nlevs);

    // Rebuild all registered MultiFabs and field aliases for pre-existing non-zero
    // levels whose BoxArrays actually changed. Guarding on `changed` makes the
    // zero-tagging case a true no-op: fill_registered_mf_level reallocates each
    // registered MF and reloads only the valid region (InterpFromCoarseLevel/
    // ParallelCopy fill no ghost cells; FillBoundary does not touch physical or
    // coarse-fine ghosts). Those halos are owned by ghostcell::start1-4, which
    // regrid does not call -- so an unconditional rebuild leaves fresh, unfilled
    // ghost bands that the next predictor/rebalance reads as garbage.
    // Level 0 is always fixed; only levels 1..min(old,new)-1 can change.
    const auto ghost_vec = amrex::IntVect(p->margin, p->margin, p->margin);
    const int rebuild_hi = std::min(old_nlevs, new_nlevs) - 1;
    if(p->mpirank==0 && changed)
        std::cout << "Rebuilding registered MultiFabs and field aliases for levels 1.." << rebuild_hi << std::endl;
    if (changed)
    for (int lev = 1; lev <= rebuild_hi; ++lev)
    {
        redefine_registered_mf_level(lev);
        rebuild_registered_field_aliases_level(lev);
        amr_cell_mf[lev].define(amrex_box_array[lev],
                                amrex_distribution_mapping[lev], 1, p->margin);
    }

    // Add new levels when the level count increased.
    if (new_nlevs > old_nlevs)
    {
        resize_registered_mf(old_nlevs, new_nlevs);
        extend_registered_fields(new_nlevs);

        amr_cell_mf.resize(new_nlevs);
        for (int lev = old_nlevs; lev < new_nlevs; ++lev)
            amr_cell_mf[lev].define(amrex_box_array[lev],
                                    amrex_distribution_mapping[lev], 1, p->margin);
    }

    // Recompute fine-mask for all active levels whenever box arrays changed.
    // Only process levels that were actually created (0 to new_nlevs-1).
    if (new_nlevs >= 2)
    {
        // Ensure amr_cell_mf is sized correctly for the new level count
        if ((int)amr_cell_mf.size() < new_nlevs)
            amr_cell_mf.resize(new_nlevs);

        amr_cell_mf[new_nlevs - 1].setVal(0);
        for (int lev = new_nlevs - 2; lev >= 0; --lev)
            amr_cell_mf[lev] = amrex::makeFineMask(
                amr_cell_mf[lev], amr_cell_mf[lev+1],
                ghost_vec, ref_vec, amrex_geometry[lev].periodicity(), 1, 0);
    }

    nlevs = new_nlevs;
    output_amrex_level_info();

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
            if (i + 1 < (gknox + 2*marge)*pow(ref_ratio,lev))
                DXP[idx] = XP[idx+1]-XP[idx];
        }
        if ((gknox + 2*marge)*pow(ref_ratio,lev) > 0)
            DXP[(gknox + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_i] = DXN[(gknox + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_i];
        for (int j = 0; j < (gknoy + 2*marge)*pow(ref_ratio,lev); ++j)
        {
            int idx = j + lev*max_j;
            DYN[idx] = YN[idx+1]-YN[idx];
            if (j + 1 < (gknoy + 2*marge)*pow(ref_ratio,lev))
                DYP[idx] = YP[idx+1]-YP[idx];
        }
        if ((gknoy + 2*marge)*pow(ref_ratio,lev) > 0)
            DYP[(gknoy + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_j] = DYN[(gknoy + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_j];
        for (int k = 0; k < (gknoz + 2*marge)*pow(ref_ratio,lev); ++k)
        {
            int idx = k + lev*max_k;
            DZN[idx] = ZN[idx+1]-ZN[idx];
            if (k + 1 < (gknoz + 2*marge)*pow(ref_ratio,lev))
                DZP[idx] = ZP[idx+1]-ZP[idx];
        }
        if ((gknoz + 2*marge)*pow(ref_ratio,lev) > 0)
            DZP[(gknoz + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_k] = DZN[(gknoz + 2*marge)*pow(ref_ratio,lev)-1 + lev*max_k];
    }
}

void grid_amrex::fill_registered_mf_level(int lev)
{
    // Ensure the BoxArray for this level is valid before using it
    if ((int)amrex_box_array.size() <= lev || amrex_box_array[lev].empty())
    {
        if(amrex::ParallelDescriptor::MyProc() == 0)
            std::cerr << "WARNING: fill_registered_mf_level called with invalid BoxArray at level "
                      << lev << std::endl;
        return;
    }

    for (auto& e : mf_registry)
    {
        if ((int)e.mf->size() <= lev) continue;

        const amrex::MultiFab& coarse_mf = (*e.mf)[lev-1];
        amrex::MultiFab&       fine_mf   = (*e.mf)[lev];

        // Preserve existing fine-level data before redefining the MultiFab.
        amrex::MultiFab old_mf;
        if (fine_mf.nComp() > 0)
            old_mf = std::move(fine_mf);

        fine_mf.define(amrex_box_array[lev],
                            amrex_distribution_mapping[lev],
                            e.ncomp, margin);

        // AMReX define()s Fabs to signalling NaN. InterpFromCoarseLevel/ParallelCopy
        // below fill only the valid region, and the trailing FillBoundary fills only
        // same-level/periodic ghosts -- the physical and coarse-fine ghost rings stay
        // NaN until ghostcell::start1-4 runs, and start1-3 do not fill a staggered
        // field's C-F face. A predictor divergence reading that NaN velocity ghost at a
        // patch low-i C-F face poisons the Poisson RHS (rhs=nan at a rest-state cell).
        // Zero the whole MF first so any ghost not reached by start1-4 is a finite 0.
        fine_mf.setVal(0.0);

        amrex::Vector<amrex::BCRec> bcrecs(e.ncomp);
        for (auto& bc : bcrecs)
            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
            {
                bc.setLo(dim, amrex::BCType::foextrap);
                bc.setHi(dim, amrex::BCType::foextrap);
            }
        amrex::PhysBCFunctNoOp null_bc;

        if (e.location >= 1 && e.location <= 3)
        {
            // Face-centred field (location 1/2/3 = FACE_X/Y/Z). The staggered velocity is held
            // in a cell-centred MultiFab with vel(cell i) == the high face == face(i+e), so a
            // cell_cons_interp would mis-place the C-F normal velocity by ~half a fine cell.
            // Convert to a face-typed temporary, prolong with face_linear_interp (the staggered,
            // divergence-aware interpolation), then copy back to cell-centred storage.
            // Assumes every component of this MF shares the same face direction.
            const int dir = e.location - 1;
            const amrex::IntVect ev = amrex::IntVect::TheDimensionVector(dir);

            amrex::MultiFab cface(amrex::convert(coarse_mf.boxArray(), ev),
                                  coarse_mf.DistributionMap(), e.ncomp, coarse_mf.nGrow());
            amrex::MultiFab fface(amrex::convert(amrex_box_array[lev], ev),
                                  amrex_distribution_mapping[lev], e.ncomp, 0);

            // cell-centred storage -> coarse face temp (face f = vel(f-e))
            for (amrex::MFIter mfi(cface); mfi.isValid(); ++mfi)
            {
                const amrex::Box& fbx = mfi.validbox();
                auto       cf = cface.array(mfi);
                const auto cc = coarse_mf.const_array(mfi);
                amrex::LoopOnCpu(fbx, e.ncomp, [&] (int i, int j, int k, int n) noexcept
                { cf(i,j,k,n) = cc(i-ev[0], j-ev[1], k-ev[2], n); });
            }

            amrex::InterpFromCoarseLevel(
                fface, 0.0,
                cface, 0, 0, e.ncomp,
                amrex_geometry[lev-1], amrex_geometry[lev],
                null_bc, 0, null_bc, 0,
                ref_vec, &amrex::face_linear_interp,
                bcrecs, 0);

            // fine face temp -> cell-centred storage (vel(i) = face(i+e))
            for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                auto       fc = fine_mf.array(mfi);
                const auto ff = fface.const_array(mfi);
                amrex::LoopOnCpu(bx, e.ncomp, [&] (int i, int j, int k, int n) noexcept
                { fc(i,j,k,n) = ff(i+ev[0], j+ev[1], k+ev[2], n); });
            }
        }
        else
        {
            amrex::InterpFromCoarseLevel(
                fine_mf, 0.0,
                coarse_mf, 0, 0, e.ncomp,
                amrex_geometry[lev-1], amrex_geometry[lev],
                null_bc, 0, null_bc, 0,
                ref_vec, &amrex::cell_cons_interp,
                bcrecs, 0);
        }

        // Overwrite interpolated cells with pre-existing fine-level data.
        // src_nghost MUST be zero: copy only from the source's VALID cells. When the
        // BoxArray is re-chopped, a destination valid cell near a moved box seam is
        // covered both by its true owner's valid cell and (within margin) by a
        // neighbouring old box's ghost cell. Feeding source ghosts (old_mf.nGrowVect())
        // makes ParallelCopy's duplicate-source resolution last-writer-wins, so valid
        // cells get clobbered by stale ghost values (extrapolated/coarse/boundary junk)
        // -> exceeding velocities on the first regrid. Valid->valid only is exact.
        if (old_mf.nComp() > 0)
        {
            fine_mf.ParallelCopy(old_mf, 0, 0, e.ncomp,
                                      amrex::IntVect::TheZeroVector(),
                                      amrex::IntVect::TheZeroVector());
            fine_mf.FillBoundary(amrex_geometry[lev].periodicity());
        }
    }

    for (auto& e : imf_registry)
    {
        if ((int)e.mf->size() <= lev) continue;

        const amrex::iMultiFab& coarse_imf = (*e.mf)[lev-1];
        amrex::iMultiFab&       fine_imf   = (*e.mf)[lev];

        // Preserve existing fine-level data before redefining the iMultiFab.
        amrex::iMultiFab old_imf;
        if (fine_imf.nComp() > 0)
            old_imf = std::move(fine_imf);

        fine_imf.define(amrex_box_array[lev],
                            amrex_distribution_mapping[lev],
                            e.ncomp, margin);

        // Zero all cells first (see the MF path): the injection + ParallelCopy below
        // fill only the valid region, so an unfilled ghost would otherwise carry an
        // undefined int flag that flag-driven loops read at the C-F ring.
        fine_imf.setVal(0);

        // Preserve the DOMAIN-EXTERIOR ghost ring (true physical walls) across the rebuild.
        // Every other cell is refilled below/after: valid cells by the coarse injection and the
        // valid->valid old copy; interior coarse-fine ghosts by flagfield's fillHigherLevels;
        // same-level/periodic ghosts by FillBoundary. But NOTHING refills the physical-boundary
        // ghost -- fillHigherLevels deliberately skips it ("must keep their OBJ_FLAG",
        // ArrayWrapper3D.cpp) assuming it survives, and the setVal(0) above just wiped it.
        // Left at 0, poisson_pcorr's wall fold/Dirichlet branches (all gated on flag4<0) never
        // fire on the fine patch's lateral domain boundary -> a stencil entry points out of the
        // patch, HYPRE drops it, the row is inconsistent and the projection injects velocity
        // (N61 blow-up). A ghost-INCLUSIVE ParallelCopy of the old data restores those exterior
        // ghosts; the seam-clobber caveat that forces src_nghost=0 on the valid copy below does
        // not apply here because every valid cell is subsequently re-set by the injection + copy.
        // (Moving-patch caveat: a patch that regrids to touch the domain boundary at a location
        //  the old patch did not will still leave those new exterior ghosts at 0; the lateral
        //  wall flag is spatially uniform, so a broadcast fallback can be added if that arises.)
        if (old_imf.nComp() > 0)
            fine_imf.ParallelCopy(old_imf, 0, 0, e.ncomp,
                                  old_imf.nGrowVect(), fine_imf.nGrowVect(),
                                  amrex_geometry[lev].periodicity());

        // Injection (piecewise-constant) from the next coarser level.
        // Build a temporary iMultiFab on the coarsened fine BoxArray so we can
        // ParallelCopy coarse data onto the fine distribution, then inject.
        const amrex::IntVect    ratio       = ref_vec;
        const int               icomp       = e.ncomp;

        amrex::BoxArray coarsened_ba = amrex::coarsen(fine_imf.boxArray(), ratio);
        amrex::iMultiFab tmp(coarsened_ba, fine_imf.DistributionMap(), icomp, 0);
        tmp.ParallelCopy(coarse_imf, 0, 0, icomp,
                         coarse_imf.nGrowVect(), amrex::IntVect::TheZeroVector());

        for (amrex::MFIter mfi(fine_imf); mfi.isValid(); ++mfi)
        {
            const amrex::Box fine_box         = mfi.validbox();
            amrex::Array4<int>       fine_arr = fine_imf.array(mfi);
            amrex::Array4<int const> tmp_arr  = tmp.const_array(mfi);
            amrex::ParallelFor(fine_box, icomp,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
                {
                    fine_arr(i, j, k, n) = tmp_arr(i/ratio[0], j/ratio[1], k/ratio[2], n);
                });
        }

        // Overwrite injected cells with pre-existing fine-level data.
        // src_nghost = 0: valid->valid only (see the MF path above -- source ghosts
        // would clobber valid cells at re-chopped box seams via ParallelCopy's
        // last-writer-wins duplicate resolution).
        if (old_imf.nComp() > 0)
        {
            fine_imf.ParallelCopy(old_imf, 0, 0, icomp,
                                  amrex::IntVect::TheZeroVector(),
                                  amrex::IntVect::TheZeroVector());
            fine_imf.FillBoundary(amrex_geometry[lev].periodicity());
        }
    }
}

void grid_amrex::resize_registered_mf(int old_nlevs, int new_nlevs)
{
    if (new_nlevs <= old_nlevs)
        return;

    for (auto& e : mf_registry)
        e.mf->resize(new_nlevs);
    for (auto& e : imf_registry)
        e.mf->resize(new_nlevs);

    for (int lev = old_nlevs; lev < new_nlevs; ++lev)
        fill_registered_mf_level(lev);
}

void grid_amrex::extend_registered_fields(int new_nlevs)
{
    for (auto* field : field_registry)
        if (field)
            field->extend_levels(new_nlevs);
}

void grid_amrex::redefine_registered_mf_level(int lev)
{
    fill_registered_mf_level(lev);
}

void grid_amrex::rebuild_registered_field_aliases_level(int lev)
{
    for (auto* field : field_registry)
        if (field)
            field->rebuild_alias_level(lev);  // no-op + cache invalidation for owning mode
}

void grid_amrex::update_registered_weno(int new_nlevs)
{
    for (auto* w : weno3_registry)
        if (w) { w->iniflag = false; }
    for (auto* w : weno5_registry)
        if (w) { w->iniflag = false; }

    lexer* p = static_cast<lexer*>(this);
    for (auto* w : weno3_registry)
        if (w) w->rebuild_levels(p, new_nlevs);
    for (auto* w : weno5_registry)
        if (w) w->rebuild_levels(p, new_nlevs);
}

void grid_amrex::create_slice_BoxArray_and_DistributionMapping(int lev)
{
    amrex::BoxList bl_slice;

    for (int i = 0; i < amrex_box_array[lev].size(); ++i)
    {
        amrex::Box b = amrex_box_array[lev][i];

        b.setSmall(2, 0);
        b.setBig(2, 0);
        bl_slice.push_back(b);
    }
    bl_slice.simplify(true);

    amrex::BoxArray ba_slice(bl_slice);
    amrex::DistributionMapping dm_slice(ba_slice);
    amrex_slice_box_array[lev] = ba_slice;
    amrex_slice_distribution_mapping[lev] = dm_slice;
}

#endif
