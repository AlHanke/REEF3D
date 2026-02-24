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

#ifndef ARRAYWRAPPER_INT_H_
#define ARRAYWRAPPER_INT_H_

#include "ArrayWrapper_base.h"

#if USE_AMREX
#include <AMReX_iMultiFab.H>
#else
#include <vector>
#endif

class lexer;

class ArrayWrapper_int : public ArrayWrapper_base<int>
{
public:
    ArrayWrapper_int(lexer* p);
    virtual ~ArrayWrapper_int() = default;

    void resize(int default_value = 0) override;

    int& operator[] (int index) override;
    operator int* () override;
    void setVal(int val, bool includeGhost = false) override;

    #if USE_AMREX
    void fillBoundary() override;
    void fillHigherLevels() override;
    #endif

private:
    #if USE_AMREX
    std::vector<amrex::iMultiFab> data;
    #else
    std::vector<std::vector<int>> data;
    #endif
    lexer* p;
};

#endif
