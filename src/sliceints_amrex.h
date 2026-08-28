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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef SLICEINTS_AMREX_H_
#define SLICEINTS_AMREX_H_

#include "sliceint_amrex.h"

class sliceint1 : public sliceint_amrex
{
public:
    sliceint1(lexer *p) : sliceint_amrex(p) {};
    virtual ~sliceint1() = default;
};

class sliceint2 : public sliceint_amrex
{
public:
    sliceint2(lexer *p) : sliceint_amrex(p) {};
    virtual ~sliceint2() = default;
};

class sliceint4 : public sliceint_amrex
{
public:
    sliceint4(lexer *p) : sliceint_amrex(p) {};
    virtual ~sliceint4() = default;
};

#endif
#endif