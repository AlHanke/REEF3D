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

#if USE_AMREX
#ifndef FIELDINT_AMREX_H_
#define FIELDINT_AMREX_H_

#include "fieldint.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

class fieldint_amrex : public fieldint
{
public:
    virtual ~fieldint_amrex() = default;

    inline int& operator()(int ii, int jj, int kk) override;

    inline int& operator()(const amrex::IntVect& iv, int comp = 0) override;

    void setVal(int val, bool includeGhost = false) override;

    void FillBoundary() override;

    amrex::iMultiFab& GetMultiFab() {return mf[p->level];};

protected:
    fieldint_amrex(lexer* p);

    lexer *p = nullptr;
    amrex::Vector<amrex::iMultiFab> mf = {};
};

#endif
#endif
