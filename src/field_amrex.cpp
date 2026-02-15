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
#include <AMReX_MultiFab.H>
#include <AMReX_Array4.H>
#include <AMReX_BCUtil.H>

field_amrex::field_amrex(lexer* p)
{
    pp = p;
}

double& field_amrex::operator()(int ii, int jj, int kk) noexcept
{
    using namespace amrex;

    IntVect cell_index{AMREX_D_DECL(ii + pp->origin_i, jj + pp->origin_j, kk + pp->origin_k)};

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        // const Box& box = mfi.fabbox();
        // FArrayBox& fab = mf[mfi];
        // Array4<Real> const& a = fab.array();
        // const auto cell_index = box.smallEnd() + IntVect{AMREX_D_DECL(ii + pp->margin, jj + pp->margin, kk + pp->margin)};
        // if(cell_index != IntVect{AMREX_D_DECL(ii + pp->origin_i, jj + pp->origin_j, kk + pp->origin_k)})
        //     Abort("field_amrex::operator(): indices don't match.");
        // const auto hi = ubound(box);
        // return a(cell_index, 0);
        if (!mfi.fabbox().contains(cell_index)) continue;
        Array4<Real> arr = mf.array(mfi);
        return arr(cell_index, 0);
    }

    Abort("field_amrex::operator(): index outside owned boxes.");
}

void field_amrex::setVal(double val, bool includeGhost)
{
    mf.setVal(val, includeGhost ? mf.nGrowVect() : amrex::IntVect{0});
}

void field_amrex::FillBoundary()
{
    mf.FillBoundary(pp->amrex_geometry.periodicity());
}

void field_amrex::FillDomainBoundary()
{
    amrex::FillDomainBoundary(mf, pp->amrex_geometry, bc);
}

void field_amrex::initialize_bc()
{
    using namespace amrex;

    bc.resize(mf.n_comp);
    for (int n = 0; n < mf.nComp(); ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            if (pp->amrex_geometry.isPeriodic(idim))
            {
                bc[n].setLo(idim, BCType::int_dir); // interior
                bc[n].setHi(idim, BCType::int_dir);
            }
            else
            {
                // ToDo: Fix this
                bc[n].setLo(idim, BCType::bogus);
                bc[n].setHi(idim, BCType::bogus);
            }
        }
    }
}
#endif
