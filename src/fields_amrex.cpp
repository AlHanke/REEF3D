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
#include "fields_amrex.h"
#include "lexer.h"
#include <AMReX_BC_TYPES.H>

// ===========================================================================
// field1
// ===========================================================================

void field1::init_params(lexer* p)
{
    auto [u_label, v_label] = field_amrex_detail::compute_parallel_u_v(p);
    params.gclabel_u    = u_label;
    params.orth_label   = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::U);
    params.awa_label    = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
    params.i10_enabled  = field_amrex_detail::compute_i10_enabled(p);
}

field1::field1(lexer* p) : field_amrex(p, DataLocation::FACE_X)
{
    init_params(p);
}

field1::field1(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp)
    : field_amrex(p, shared_mf, comp, DataLocation::FACE_X)
{
    init_params(p);
}

void field1::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field1BcDecision(params));
}

void field1::UpdateBCRecs(int gcv)
{
    UpdateBCRecsImpl(gcv, amrex_bc_func::Field1BcDecision(params));
}

// ===========================================================================
// field2
// ===========================================================================

void field2::init_params(lexer* p)
{
    auto [u_label, v_label] = field_amrex_detail::compute_parallel_u_v(p);
    params.gclabel_v    = v_label;
    params.orth_label   = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::V);
    params.awa_label    = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
}

field2::field2(lexer* p) : field_amrex(p, DataLocation::FACE_Y)
{
    init_params(p);
}

field2::field2(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp)
    : field_amrex(p, shared_mf, comp, DataLocation::FACE_Y)
{
    init_params(p);
}

void field2::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field2BcDecision(params));
}

void field2::UpdateBCRecs(int gcv)
{
    UpdateBCRecsImpl(gcv, amrex_bc_func::Field2BcDecision(params));
}

// ===========================================================================
// field3
// ===========================================================================

void field3::init_params(lexer* p)
{
    params.gclabel_w    = field_amrex_detail::compute_parallel_w(p);
    params.orth_label   = field_amrex_detail::compute_orth_label(p);
    params.inflow_label = field_amrex_detail::compute_inflow_label(p);
    params.outflow_label = field_amrex_detail::compute_outflow_label(p, field_amrex_detail::OutflowAxis::W);
    params.awa_label    = field_amrex_detail::compute_awa_label(p);
    params.gclabel_outflow = field_amrex_detail::compute_gclabel_outflow(p);
    params.A10 = p->A10;
}

field3::field3(lexer* p) : field_amrex(p, DataLocation::FACE_Z)
{
    init_params(p);
}

field3::field3(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp)
    : field_amrex(p, shared_mf, comp, DataLocation::FACE_Z)
{
    init_params(p);
}

void field3::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field3BcDecision(params));
}

void field3::UpdateBCRecs(int gcv)
{
    UpdateBCRecsImpl(gcv, amrex_bc_func::Field3BcDecision(params));
}

// ===========================================================================
// field4
// ===========================================================================

// Shared by field4 and field7: Field7BcDecision aliases Field4BcDecision, so
// both fill the same Field4Params. Factored out rather than duplicated so the
// two cannot drift apart.
static void init_field4_params(lexer* p,
                               amrex_bc_func::Field4BcDecision::Field4Params& params)
{
    params.B77  = p->B77;
    params.H61  = p->H61;
    params.H62  = p->H62;
    params.H63  = p->H63;
    params.H64  = p->H64;
    params.H65  = p->H65;
    params.H66  = p->H66;
    params.pressout_label = (p->B77 == 1 || p->B77 == 10);
    params.pressin_label  = (p->B76 != 1);
    params.awa_label = field_amrex_detail::compute_awa_label(p);
    params.gclabel_lsm_in_neumann  = !(p->I230 > 0 || p->B98 >= 3 || p->B60 > 0);
    params.gclabel_press_in_neumann = (p->B76 != 2 && p->B76 != 3);
}

void field4::init_params(lexer* p)
{
    init_field4_params(p, params);
}

field4::field4(lexer* p) : field_amrex(p, DataLocation::CELL_CENTERED)
{
    init_params(p);
}

field4::field4(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp)
    : field_amrex(p, shared_mf, comp, DataLocation::CELL_CENTERED)
{
    init_params(p);
}

void field4::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field4BcDecision(params));
}

void field4::UpdateBCRecs(int gcv)
{
    UpdateBCRecsImpl(gcv, amrex_bc_func::Field4BcDecision(params));
}

// ===========================================================================
// field7 — sigma-grid vertical nodes (z-nodal MultiFab)
// ===========================================================================

void field7::init_params(lexer* p)
{
    init_field4_params(p, params);
}

field7::field7(lexer* p) : field_amrex(p, DataLocation::NODE_Z)
{
    init_params(p);
}

field7::field7(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp)
    : field_amrex(p, shared_mf, comp, DataLocation::NODE_Z)
{
    init_params(p);
}

void field7::FillDomainBoundary(int gcv)
{
    FillDomainBoundaryImpl(gcv, amrex_bc_func::Field7BcDecision(params));
}

void field7::UpdateBCRecs(int gcv)
{
    UpdateBCRecsImpl(gcv, amrex_bc_func::Field7BcDecision(params));
}
#endif
