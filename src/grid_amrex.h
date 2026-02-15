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

    // AMReX Geometry
    std::vector<amrex::Geometry> amrex_geometry;
    std::vector<amrex::BoxArray> amrex_box_array;
    std::vector<amrex::DistributionMapping> amrex_distribution_mapping;
    std::vector<amrex::MultiFab> amr_mf;
    std::unique_ptr<amrex::MFIter> default_mfi;
    amrex::MFIter* amr_mfi = nullptr;
    std::vector<amrex::iMultiFab> flag1_imf;
    std::vector<amrex::iMultiFab> flag2_imf;
    std::vector<amrex::iMultiFab> flag3_imf;
    std::vector<amrex::iMultiFab> flag4_imf;
    std::vector<amrex::iMultiFab> flag7_imf;

    int level;
    const int nlevs = 1;
    const int ncomp = 1;
    int bc_type[6] = {0,0,0,0,0,0};

protected:
    void setup_amrex_geometry(lexer* p, ghostcell* pgc);
};

#endif
#endif
