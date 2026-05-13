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
#ifndef GRID_AMREX_H_
#define GRID_AMREX_H_

#include "grid.h"

// AMReX includes
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

#include <memory>
#include <vector>
#include <utility>
#include <algorithm>

class lexer;
class ghostcell;
class fdm;

namespace amrex {
    class Geometry;
    class BoxArray;
    class DistributionMapping;
}

/*!
    * @brief grid_amrex is a class that incorporates AMReX data structures and functionalities for handling the computational grid in REEF3D.
    * It inherits from the grid class and extends it with AMReX-specific members and methods.
    * The class includes methods for defining inflow and outflow boundary areas, setting up AMReX geometry, and managing AMReX MultiFabs and iMultiFabs for cell-centered data and flags.
    * It also includes members for storing AMReX Geometry, BoxArray, and DistributionMapping for each level of the AMR hierarchy.
    * The class is designed to facilitate the integration of AMReX's capabilities for handling complex grid structures and parallel computations within the REEF3D framework.
*/
class grid_amrex : public grid
{
public:
    grid_amrex() = default;
    virtual ~grid_amrex() = default;

    void define_inflow_outflow_ba();
    void regrid_amrex_box_array_and_distribution_mapping(lexer* p, fdm* a);

    void update_cell_coordinates();

    // Registry for all owning MultiFab/iMultiFab vectors so they can be
    // resized collectively during AMR regrid.
    struct MFEntry  { amrex::Vector<amrex::MultiFab>*  mf; int ncomp; };
    struct IMFEntry { amrex::Vector<amrex::iMultiFab>* mf; int ncomp; };
    std::vector<MFEntry>  mf_registry;
    std::vector<IMFEntry> imf_registry;

    void register_mf (amrex::Vector<amrex::MultiFab>*  mf, int ncomp) { mf_registry .push_back({mf, ncomp}); }
    void register_imf(amrex::Vector<amrex::iMultiFab>* mf, int ncomp) { imf_registry.push_back({mf, ncomp});}

    void deregister_mf(amrex::Vector<amrex::MultiFab>* mf)
    {
        auto it = std::remove_if(mf_registry.begin(), mf_registry.end(),
                                 [mf](const MFEntry& entry) { return entry.mf == mf; });
        mf_registry.erase(it, mf_registry.end());
    }
    void deregister_imf(amrex::Vector<amrex::iMultiFab>* mf)
    {
        auto it = std::remove_if(imf_registry.begin(), imf_registry.end(),
                                 [mf](const IMFEntry& entry) { return entry.mf == mf; });
        imf_registry.erase(it, imf_registry.end());
    }

    // Registry for field_amrex objects so they can extend BCRecs/m_alias on regrid
    std::vector<class field_amrex*> field_registry;
    void register_field(class field_amrex* f) { field_registry.push_back(f); }
    void extend_registered_fields(int new_nlevs);
    void deregister_field(class field_amrex* f)
    {
        auto it = std::remove(field_registry.begin(), field_registry.end(), f);
        field_registry.erase(it, field_registry.end());
    }

    // Registry for WENO non-uniform-grid objects so their precomputed tables
    // are reallocated and recalculated alongside coordinate arrays on regrid.
    std::vector<class weno3_nug_func*> weno3_registry;
    std::vector<class weno_nug_func*>  weno5_registry;

    void register_weno3(class weno3_nug_func* w) { weno3_registry.push_back(w); }
    void register_weno5(class weno_nug_func*  w) { weno5_registry.push_back(w); }
    void deregister_weno3(class weno3_nug_func* w)
    {
        auto it = std::remove(weno3_registry.begin(), weno3_registry.end(), w);
        weno3_registry.erase(it, weno3_registry.end());
    }
    void deregister_weno5(class weno_nug_func* w)
    {
        auto it = std::remove(weno5_registry.begin(), weno5_registry.end(), w);
        weno5_registry.erase(it, weno5_registry.end());
    }
    void update_registered_weno(int new_nlevs);

    void resize_registered_mf(int old_nlevs, int new_nlevs);
    void fill_registered_mf_level(int lev);
    /// Redefine all registered MultiFabs at @p lev with the current BoxArray and
    /// DistributionMapping, filling data via coarse-level interpolation.
    /// Called for existing non-zero levels whose box arrays change during regrid.
    void redefine_registered_mf_level(int lev);
    /// Rebuild view-mode aliases at @p lev for all registered field_amrex objects.
    void rebuild_registered_field_aliases_level(int lev);
    void update_cell_spacing();

    // AMReX Data structures
    amrex::Vector<amrex::Geometry> amrex_geometry; // Phyiscal domain and coordinate system
    amrex::Vector<amrex::BoxArray> amrex_box_array; // BoxArray defines the index space decomposition of the domain into boxes
    amrex::Vector<amrex::DistributionMapping> amrex_distribution_mapping; // DistributionMapping defines the mapping of boxes in the BoxArray to MPI ranks
    amrex::Vector<amrex::Vector<std::pair<amrex::RealVect,amrex::RealVect>>> amrex_refined_grid_coords; // Input: Coordinates of the refined grid boxes for each level, index is offset by 1 (i.e. amrex_refined_grid_coords[0] is for level 1, etc.)

    // Looping structures
    amrex::Vector<amrex::iMultiFab> amr_cell_mf;
    std::unique_ptr<amrex::MFIter> default_cell_mfi;
    amrex::MFIter* amr_cell_mfi = nullptr;

    // Tile bounds — set once per tile by set_tile_mfi; shared with field_amrex cache.
    amrex::Dim3 amr_tile_lo    = {0, 0, 0};
    amrex::Dim3 amr_tile_hi    = {0, 0, 0};
    int amr_fab_mfi_idx        = -1; // Index of current MFIter
    int amr_local_tile_idx     = -1; // Local tile index within the current MFIter

    /// Set the active MFIter, compute and cache the tile's lo/hi bounds and FAB
    /// index, then return the previously active MFIter (for guard-struct restore).
    /// Called exactly once per tile loop iteration by TILE_LOOP.
    inline amrex::MFIter* set_tile_mfi(amrex::MFIter* mfi) noexcept
    {
        amrex::MFIter* old = amr_cell_mfi;
        amr_cell_mfi = mfi;

        const amrex::Box tb = mfi->tilebox();
        amr_tile_lo         = amrex::lbound(tb);
        amr_tile_hi         = amrex::ubound(tb);
        amr_fab_mfi_idx     = mfi->index();
        amr_local_tile_idx  = mfi->LocalTileIndex();

        return old;
    }

    // Inflow and outflow areas
    amrex::Vector<amrex::iMultiFab> inflow_ba;
    amrex::Vector<amrex::Vector<amrex::IntVect>> inflow_ijk;
    amrex::Vector<amrex::iMultiFab> outflow_ba;
    amrex::Vector<amrex::Vector<amrex::IntVect>> outflow_ijk;

    int nlevs = 1; // Number of AMR levels
    amrex::IntVect ref_vec;
    static constexpr int ref_ratio = 2;
    const int ncomp = 1;
    int bc_type[6] = {0,0,0,0,0,0};

protected:
    void setup_amrex_geometry(lexer* p, ghostcell* pgc);

    static constexpr int max_nlevs = 3;
private:
    void create_amrex_box_array_and_distribution_mapping_level_n();
    void output_amrex_level_info();
};

#endif
#endif
