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

#include "fieldints_amrex.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>

fieldint1::fieldint1(lexer *p) : fieldint_amrex(p)
{
    LevelLOOP
    {
        mf[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, p->margin);
    }
}

fieldint2::fieldint2(lexer *p) : fieldint_amrex(p)
{
    LevelLOOP
    {
        mf[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, p->margin);
    }
}

fieldint3::fieldint3(lexer *p) : fieldint_amrex(p)
{
    LevelLOOP
    {
        mf[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, p->margin);
    }
}

fieldint4::fieldint4(lexer *p) : fieldint_amrex(p)
{
    LevelLOOP
    {
        mf[p->level].define(p->amrex_box_array[p->level], p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, p->margin);
    }
}
