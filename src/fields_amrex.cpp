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
Author: Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#if USE_AMREX
#include "fields_amrex.h"
#include "lexer.h"
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BC_TYPES.H>

field1::field1(lexer *p) : field_amrex(p, amrex_bc_func::DataLocation::FACE_X)
{
    LevelLOOP
    {
        amrex::BoxArray box = p->amrex_box_array[p->level];
        // box = amrex::convert(p->amrex_box_array[p->level], amrex::IntVect{AMREX_D_DECL(1,0,0)});
        // box.grow(amrex::IntVect{AMREX_D_DECL(-1,0,0)});
        mf[p->level].define(box, p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, mf[p->level].nGrow());
    }

    auto [u_label, v_label] = field_amrex_detail::compute_parallel_u_v(p);
    params.gclabel_u = u_label;
    params.orth_label = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::U);
    params.awa_label = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
    params.i10_enabled = field_amrex_detail::compute_i10_enabled(p);
}

field2::field2(lexer *p) : field_amrex(p, amrex_bc_func::DataLocation::FACE_Y)
{
    LevelLOOP
    {
        amrex::BoxArray box = p->amrex_box_array[p->level];
        // box = amrex::convert(p->amrex_box_array[p->level], amrex::IntVect{AMREX_D_DECL(0,1,0)});
        // box.grow(amrex::IntVect{AMREX_D_DECL(0,-1,0)});
        mf[p->level].define(box, p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, mf[p->level].nGrow());
    }

    auto [u_label, v_label] = field_amrex_detail::compute_parallel_u_v(p);
    params.gclabel_v = v_label;
    params.orth_label = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::V);
    params.awa_label = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
}

field3::field3(lexer *p) : field_amrex(p, amrex_bc_func::DataLocation::FACE_Z)
{
    LevelLOOP
    {
        amrex::BoxArray box = p->amrex_box_array[p->level];
        // box = amrex::convert(p->amrex_box_array[p->level], amrex::IntVect{AMREX_D_DECL(0,0,1)});
        // box.grow(amrex::IntVect{AMREX_D_DECL(0,0,-1)});
        mf[p->level].define(box, p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, mf[p->level].nGrow());
    }

    params.gclabel_w = field_amrex_detail::compute_parallel_w(p);
    params.orth_label = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::W);
    params.awa_label = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
    params.A10 = p->A10;
}

field4::field4(lexer *p) : field_amrex(p, amrex_bc_func::DataLocation::CELL_CENTERED)
{
    LevelLOOP
    {
        amrex::BoxArray box = p->amrex_box_array[p->level];
        mf[p->level].define(box, p->amrex_distribution_mapping[p->level], p->ncomp, p->margin);
        mf[p->level].setVal(0, 0, mf[p->level].n_comp, mf[p->level].nGrow());
    }

    params.B77 = p->B77;
    params.H61 = p->H61;
    params.H62 = p->H62;
    params.H63 = p->H63;
    params.H64 = p->H64;
    params.H65 = p->H65;
    params.H66 = p->H66;
    params.pressout_label = (p->B77 == 1 || p->B77 == 10);
    params.pressin_label = (p->B76 != 1);
    params.awa_label = field_amrex_detail::compute_awa_label(p);
    params.gclabel_lsm_in_neumann = !(p->I230 > 0 || p->B98 >= 3 || p->B60 > 0);
    params.gclabel_press_in_neumann = (p->B76 != 2 && p->B76 != 3);
}

void field1::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field1BcDecision(params));
}

void field2::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field2BcDecision(params));
}

void field3::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field3BcDecision(params));
}

void field4::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field4BcDecision(params));
}

#endif
