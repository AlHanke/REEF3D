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
#ifndef FIELDINTS_AMReX_H_
#define FIELDINTS_AMReX_H_

#include"fieldint_amrex.h"

class fieldint1 : public fieldint_amrex
{
public:
    fieldint1 (lexer*);
    virtual ~fieldint1() = default;
};

class fieldint2 : public fieldint_amrex
{
public:
    fieldint2 (lexer*);
    virtual ~fieldint2() = default;
};

class fieldint3 : public fieldint_amrex
{
public:
    fieldint3 (lexer*);
    virtual ~fieldint3() = default;
};

class fieldint4 : public fieldint_amrex
{
public:
    fieldint4 (lexer*);
    virtual ~fieldint4() = default;
};

#endif
#endif
