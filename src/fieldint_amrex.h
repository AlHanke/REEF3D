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
#ifndef FIELDINT_AMREX_H_
#define FIELDINT_AMREX_H_

#include "fieldint.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>

class fieldint_amrex : public fieldint
{
public:
    virtual ~fieldint_amrex() = default;

    int& operator()(int, int, int) noexcept override final;

    void FillBoundary();

protected:
    fieldint_amrex(lexer* p);

    lexer *pp;
    amrex::iMultiFab mf;
    const int num_components = 1;
};

#endif
#endif
