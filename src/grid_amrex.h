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

class grid_amrex : public grid
{
public:
    grid_amrex() = default;
    virtual ~grid_amrex() = default;

    void define_inflow_outflow_ba();

    // AMReX Geometry
    amrex::Vector<amrex::Geometry> amrex_geometry;
    amrex::Vector<amrex::BoxArray> amrex_box_array;
    amrex::Vector<amrex::DistributionMapping> amrex_distribution_mapping;

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

    int level;
    const int nlevs = 1;
    const int ncomp = 1;
    int bc_type[6] = {0,0,0,0,0,0};

protected:
    void setup_amrex_geometry(lexer* p, ghostcell* pgc);
};

#endif
