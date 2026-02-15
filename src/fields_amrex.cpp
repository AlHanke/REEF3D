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

#include"fields_amrex.h"
#include"lexer.h"
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>

field1_amrex::field1_amrex(lexer *p) : field_amrex(p)
{
    amrex::BoxArray box = p->amrex_box_array;
    // box = amrex::convert(p->amrex_box_array, amrex::IntVect{AMREX_D_DECL(1,0,0)});
    mf.define(box, p->amrex_distribution_mapping, num_components, p->margin);
    mf.setVal(0, 0, mf.n_comp, p->margin);
    mf.setVal(0);
    fillBoundary();
}

field2_amrex::field2_amrex(lexer *p) : field_amrex(p)
{
    amrex::BoxArray box = p->amrex_box_array;
    // box = amrex::convert(p->amrex_box_array, amrex::IntVect{AMREX_D_DECL(0,1,0)});
    mf.define(box, p->amrex_distribution_mapping, num_components, p->margin);
    mf.setVal(0, 0, mf.n_comp, p->margin);
    mf.setVal(0);
    fillBoundary();
}

field3_amrex::field3_amrex(lexer *p) : field_amrex(p)
{
    amrex::BoxArray box = p->amrex_box_array;
    // box = amrex::convert(p->amrex_box_array, amrex::IntVect{AMREX_D_DECL(0,0,1)});
    mf.define(box, p->amrex_distribution_mapping, num_components, p->margin);
    mf.setVal(0, 0, mf.n_comp, p->margin);
    mf.setVal(0);
    fillBoundary();
}

field4_amrex::field4_amrex(lexer *p) : field_amrex(p)
{
    amrex::BoxArray box = p->amrex_box_array;
    // box = amrex::convert(p->amrex_box_array, amrex::IntVect{AMREX_D_DECL(0,0,0)});
    mf.define(box, p->amrex_distribution_mapping, num_components, p->margin);
    mf.setVal(0, 0, mf.n_comp, p->margin);
    mf.setVal(0);
    fillBoundary();
}
