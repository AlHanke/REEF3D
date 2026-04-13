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

#include "fieldint_amrex.h"
#include "lexer.h"
#include <AMReX_BCUtil.H>
#include <AMReX_MFIter.H>
#include <AMReX_iMultiFab.H>

fieldint_amrex::fieldint_amrex(lexer* p)
{
    pp = p;
}

int& fieldint_amrex::operator()(int ii, int jj, int kk)
{
    using namespace amrex;

    IntVect cell_index{AMREX_D_DECL(ii + pp->origin_i, jj + pp->origin_j, kk + pp->origin_k)};

    for (MFIter mfi(mf[pp->level]); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.fabbox();
        const auto cell_index = box.smallEnd() + IntVect{AMREX_D_DECL(ii + pp->margin, jj + pp->margin, kk + pp->margin)};
        return mf[pp->level].array(mfi)(cell_index, 0);
    }
    // return (mf[pp->level][*(pp->amr_mfi)].array()(amrex::IntVect{AMREX_D_DECL(ii, jj, kk)} + amrex::IntVect{amrex::lbound(pp->amr_mfi->validbox())}, 0));
}

void fieldint_amrex::setVal(int val, bool includeGhost)
{
    mf[pp->level].setVal(val, (includeGhost ? pp->margin : 0));
}

void fieldint_amrex::fillBoundary()
{
    mf[pp->level].FillBoundary(pp->amrex_geometry[pp->level].periodicity());
}

void fieldint_amrex::FillDomainBoundary()
{
    // amrex::FillDomainBoundary(mf[pp->level], pp->amrex_geometry[pp->level], bc[pp->level]);
}

void fieldint_amrex::initialize_bc()
{
    using namespace amrex;

    bc.resize(pp->nlevs);
    for(pp->level=0; pp->level<pp->nlevs; ++pp->level)
    {
        bc[pp->level].resize(mf[pp->level].n_comp);
        for (int n = 0; n < mf[pp->level].nComp(); ++n)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (pp->amrex_geometry[pp->level].isPeriodic(idim))
                {
                    bc[pp->level][n].setLo(idim, BCType::int_dir); // interior
                    bc[pp->level][n].setHi(idim, BCType::int_dir);
                }
                else
                {
                    // ToDo: Fix this
                    bc[pp->level][n].setLo(idim, BCType::bogus);
                    bc[pp->level][n].setHi(idim, BCType::bogus);
                }
            }
        }
    }
}
