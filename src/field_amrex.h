/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef FIELD_AMREX_H_
#define FIELD_AMREX_H_

#include "field.h"
#include "lexer.h"
#include <AMReX_MultiFab.H>
#include <vector>

class field_amrex : public field
{
public:
    virtual ~field_amrex() = default;

    double& operator()(int ii, int jj, int kk) noexcept override final;

    void setVal(double val, bool includeGhost = false) override final;

    void FillBoundary() override;

    void FillDomainBoundary() override;

protected:
    field_amrex(lexer* p);

    lexer *p;
    std::vector<amrex::MultiFab> mf;
};

#endif
#endif
