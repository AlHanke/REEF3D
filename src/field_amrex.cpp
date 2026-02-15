/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include "field_amrex.h"
#include "lexer.h"
#include "amrex_bc_func.h"
#include <AMReX_BCUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>

field_amrex::field_amrex(lexer* p): face_bc_values{p->bcside1, p->bcside4, p->bcside3, p->bcside2, p->bcside5, p->bcside6}
{
    field_amrex::p = p;
    mf.resize(p->nlevs);

    BCRecs.resize(p->nlevs);
    for (auto& bc_rec : BCRecs)
        bc_rec.resize(p->ncomp);
}

double& field_amrex::operator()(int ii, int jj, int kk) noexcept
{
    return (mf[p->level][*(p->amr_mfi)].array()(amrex::IntVect{AMREX_D_DECL(ii, jj, kk)} + amrex::IntVect{amrex::lbound(p->amr_mfi->validbox())}, 0));
}

void field_amrex::setVal(double val, bool includeGhost)
{
    mf[p->level].setVal(val, includeGhost ? mf[p->level].nGrowVect() : amrex::IntVect{0});
}

void field_amrex::FillBoundary()
{
    mf[p->level].FillBoundary(p->amrex_geometry[p->level].periodicity());
}

void field_amrex::FillDomainBoundary()
{
    amrex::Vector<amrex::Box> loc_boxes;
    for (amrex::MFIter mfi(mf[p->level]); mfi.isValid(); ++mfi) {
        loc_boxes.push_back(mfi.validbox());
    }

    amrex::Gpu::DeviceVector<amrex::Box> device_boxes(loc_boxes.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, loc_boxes.begin(), loc_boxes.end(), device_boxes.begin());

    amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFill> bf(amrex_bc_func::MyExtBCFill{face_bc_values, device_boxes.data(), static_cast<int>(device_boxes.size())});

    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFill>> physbcf(p->amrex_geometry[p->level], BCRecs[p->level], bf);
    amrex::FillPatchSingleLevel(mf[p->level], amrex::Real(p->simtime),
                                {&(mf[p->level])}, {amrex::Real(p->simtime)},
                                0, 0, mf[p->level].nComp(), p->amrex_geometry[p->level], physbcf, 0);
}
#endif
