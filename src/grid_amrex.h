/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke (@AlHanke)
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

class lexer;
class ghostcell;

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

    void update_cell_coordinates();
    void update_cell_spacing();

    // AMReX Data structures
    amrex::Vector<amrex::Geometry> amrex_geometry; // Phyiscal domain and coordinate system
    amrex::Vector<amrex::BoxArray> amrex_box_array; // BoxArray defines the index space decomposition of the domain into boxes
    amrex::Vector<amrex::DistributionMapping> amrex_distribution_mapping; // DistributionMapping defines the mapping of boxes in the BoxArray to MPI ranks
    amrex::Vector<amrex::Vector<std::pair<amrex::RealVect,amrex::RealVect>>> amrex_refined_grid_coords; // Input: Coordinates of the refined grid boxes for each level, index is offset by 1 (i.e. amrex_refined_grid_coords[0] is for level 1, etc.)

    // Looping structures
    amrex::Vector<amrex::MultiFab> amr_cell_mf;
    std::unique_ptr<amrex::MFIter> default_cell_mfi;
    amrex::MFIter* amr_cell_mfi = nullptr;

    // Inflow and outflow areas
    amrex::Vector<amrex::iMultiFab> inflow_ba;
    amrex::Vector<amrex::Vector<amrex::IntVect>> inflow_ijk;
    amrex::Vector<amrex::iMultiFab> outflow_ba;
    amrex::Vector<amrex::Vector<amrex::IntVect>> outflow_ijk;

    // Flags using iMultiFab
    amrex::Vector<amrex::iMultiFab> flag1_iMF;
    amrex::Vector<amrex::iMultiFab> flag2_iMF;
    amrex::Vector<amrex::iMultiFab> flag3_iMF;
    amrex::Vector<amrex::iMultiFab> flag4_iMF;
    amrex::Vector<amrex::iMultiFab> flag7_iMF;

    int nlevs = 1; // Number of AMR levels
    const int ref_ratio = 2;
    const int ncomp = 1;
    int bc_type[6] = {0,0,0,0,0,0};

protected:
    void setup_amrex_geometry(lexer* p, ghostcell* pgc);
private:
    void create_amrex_box_array_and_distribution_mapping_level_n();
};

#endif
#endif
