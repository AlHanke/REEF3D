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

#ifndef FIELD_AMREX_H_
#define FIELD_AMREX_H_

#include "field.h"
#include "lexer.h"
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

class field_amrex : public field
{
public:
    virtual ~field_amrex() = default;

    double& operator()(int, int, int) override;

    void fillBoundary() override;

    void FillDomainBoundary() override;

    amrex::MultiFab& GetMultiFab() override { return mf; }

    void setVal(double val) override { mf.setVal(val,0); }

protected:
    field_amrex(lexer* p);

    void initialize_bc();

    lexer *pp;
    amrex::MultiFab mf;
    amrex::Vector<amrex::BCRec> bc;
    const int num_components = 1;
};

#endif
