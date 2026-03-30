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

    void resize(int default_value = 0) override final;

    inline int& operator[] (int index) override final;
    operator int* () override final;
    void setVal(int val, bool includeGhost = false) override final;

    #if USE_AMREX
    void fillBoundary() override final;
    void fillHigherLevels() override final;
    inline amrex::iMultiFab& GetMultiFab();
    inline const amrex::iMultiFab& GetMultiFab() const;
    inline amrex::iMultiFab& GetMultiFab(int level) {return data[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return data[level];};
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
