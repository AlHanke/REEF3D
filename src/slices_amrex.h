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
#ifndef SLICES_AMReX_H_
#define SLICES_AMReX_H_

#include "slice_amrex.h"

class slice1_amrex : public slice_amrex
{
public:
    slice1_amrex(lexer *pp) : slice_amrex(pp, DataLocation::FACE_X) {};
    virtual ~slice1_amrex() = default;

    void FillDomainBoundary(int gcv) override final
    {
        FillDomainBoundaryImpl(gcv, amrex_bc_func2D::Slice1BcDecision(params));
    }
    void UpdateBCRecs(int gcv) override final
    {
        UpdateBCRecsImpl(gcv, amrex_bc_func2D::Slice1BcDecision(params));
    };
private:
    void init_params(lexer *p)
    {
        params.gclabel_u = p->B20==1 ? amrex_bc_func2D::BoundaryConditionTypeLabel::NEUMANN : amrex_bc_func2D::BoundaryConditionTypeLabel::NOSLIP;
        params.awa_label = p->B99 >= 3;
        params.B99 = p->B99;
    };
    amrex_bc_func2D::Slice1BcDecision::Slice1Params params{};
};

class slice2_amrex : public slice_amrex
{
public:
    slice2_amrex(lexer *pp) : slice_amrex(pp, DataLocation::FACE_Y) {};
    virtual ~slice2_amrex() = default;

    void FillDomainBoundary(int gcv) override final
    {
        FillDomainBoundaryImpl(gcv, amrex_bc_func2D::Slice2BcDecision(params));
    }
    void UpdateBCRecs(int gcv) override final
    {
        UpdateBCRecsImpl(gcv, amrex_bc_func2D::Slice2BcDecision(params));
    };
private:
    void init_params(lexer *p)
    {
        params.gclabel_v = p->B20==1 ? amrex_bc_func2D::BoundaryConditionTypeLabel::NEUMANN : amrex_bc_func2D::BoundaryConditionTypeLabel::NOSLIP;
        params.B99 = p->B99;
    };
    amrex_bc_func2D::Slice2BcDecision::Slice2Params params{};
};

class slice4_amrex : public slice_amrex
{
public:
    slice4_amrex(lexer *pp) : slice_amrex(pp, DataLocation::CELL_CENTERED) {};
    virtual ~slice4_amrex() = default;

    void FillDomainBoundary(int gcv) override final
    {
        FillDomainBoundaryImpl(gcv, amrex_bc_func2D::Slice4BcDecision(params));
    }
    void UpdateBCRecs(int gcv) override final
    {
        UpdateBCRecsImpl(gcv, amrex_bc_func2D::Slice4BcDecision(params));
    };
private:
    void init_params(lexer *p)
    {
        params.A515 = p->A515;
        params.B98 = p->B98;
        params.B99 = p->B99;
    };
    amrex_bc_func2D::Slice4BcDecision::Slice4Params params{};
};

#endif
#endif