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

field_amrex::field_amrex(lexer* p, amrex_bc_func::DataLocation data_location): const_params{{p->bcside1, p->bcside4, p->bcside3, p->bcside2, p->bcside5, p->bcside6},
                                           {p->H61_T, p->H64_T, p->H63_T, p->H62_T, p->H65_T, p->H66_T}, bool(p->j_dir), data_location}
{
    field_amrex::p = p;
    mf.resize(p->nlevs);

    BCRecs.resize(p->nlevs);
    for (auto& bc_rec : BCRecs)
        bc_rec.resize(p->ncomp);
}

double& field_amrex::operator()(int ii, int jj, int kk)
{
    return (mf[p->level][*(p->amr_cell_mfi)].array()(amrex::IntVect(AMREX_D_DECL(ii, jj, kk)) + amrex::IntVect(amrex::lbound(p->amr_cell_mfi->tilebox())), 0));
}

double& field_amrex::operator()(const amrex::IntVect& iv, int comp)
{
    return (mf[p->level][*(p->amr_cell_mfi)].array()(iv, comp));
}

void field_amrex::setVal(double val, bool includeGhost)
{
    mf[p->level].setVal(val, (includeGhost ? p->margin : 0));
}

void field_amrex::FillBoundary()
{
    mf[p->level].FillBoundary(p->amrex_geometry[p->level].periodicity());
}

void field_amrex::FillDomainBoundaryValue(double value, int dir, bool high)
{
    LevelLOOP
    {
        amrex::Box dom = p->amrex_geometry[p->level].Domain();
        TileLOOP
        {
            const amrex::Box& validbox = p->amr_cell_mfi->validbox();
            amrex::Box gbx = validbox;
            bool apply = false;

            amrex::IntVect ng(p->margin);
            if (!const_params.y_dimension_exists) ng[1] = 0;
            gbx.grow(ng);

            // Apply boundary condition at x_min
            if ((validbox.smallEnd(0) == dom.smallEnd(0)) && (dir == 0 && !high))
            {
                gbx.setBig(0, dom.smallEnd(0) - 1);
                apply = true;
            }
            // Apply boundary condition at x_max
            if ((validbox.bigEnd(0) == dom.bigEnd(0)) && (dir == 0 && high))
            {
                gbx.setSmall(0, dom.bigEnd(0) + 1);
                apply = true;
            }

            // Apply boundary condition at y_min
            if ((validbox.smallEnd(1) == dom.smallEnd(1)) && (dir == 1 && !high))
            {
                gbx.setBig(1, dom.smallEnd(1) - 1);
                apply = true;
            }
            // Apply boundary condition at y_max
            if ((validbox.bigEnd(1) == dom.bigEnd(1)) && (dir == 1 && high))
            {
                gbx.setSmall(1, dom.bigEnd(1) + 1);
                apply = true;
            }

            // Apply boundary condition at z_min
            if ((validbox.smallEnd(2) == dom.smallEnd(2)) && (dir == 2 && !high))
            {
                gbx.setBig(2, dom.smallEnd(2) - 1);
                apply = true;
            }
            // Apply boundary condition at z_max
            if ((validbox.bigEnd(2) == dom.bigEnd(2)) && (dir == 2 && high))
            {
                gbx.setSmall(2, dom.bigEnd(2) + 1);
                apply = true;
            }

            if (apply)
            {
                auto arr = mf[p->level][*(p->amr_cell_mfi)].array();
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    arr(i, j, k) = value;
                });
            }
        }
    }
}
