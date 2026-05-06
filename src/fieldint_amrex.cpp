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
#include "fieldint_amrex.h"
#include "lexer.h"
#include <AMReX_BCUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_DistributionMapping.H>

fieldint_amrex::fieldint_amrex(lexer* p)
{
    fieldint_amrex::p = p;
    mf = make_imf(p, p->ncomp);
}

void fieldint_amrex::setVal(int val, bool includeGhost)
{
    LEVEL_LOOP
    {
        mf[p->level].setVal(val, includeGhost ? mf[p->level].nGrowVect() : amrex::IntVect{0});
    }
}

void fieldint_amrex::FillBoundary()
{
    LEVEL_LOOP
    {
        mf[p->level].FillBoundary(p->amrex_geometry[p->level].periodicity());
    }
}
#endif
