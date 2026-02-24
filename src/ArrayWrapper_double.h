/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#ifndef ARRAYWRAPPER_DOUBLE_H_
#define ARRAYWRAPPER_DOUBLE_H_

#include "ArrayWrapper_base.h"
#include "lexer.h"

#if USE_AMREX
#include <AMReX_MultiFab.H>
#else
#include <vector>
#endif

class lexer;

class ArrayWrapper_double : public ArrayWrapper_base<double>
{
public:
    ArrayWrapper_double(lexer* p);
    virtual ~ArrayWrapper_double() = default;

    void resize(double default_value = 0.0) override;

    double& operator[] (int index) override;
    operator double* () override;
    void setVal(double val, bool includeGhost = false) override;

    #if USE_AMREX
    void fillBoundary() override;
    void fillHigherLevels() override;
    #endif

private:
    #if USE_AMREX
    std::vector<amrex::MultiFab> data;
    #else
    std::vector<std::vector<double>> data;
    #endif
    lexer* p;
};

#endif
