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
Author: Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#if (not USE_AMREX)
#ifndef FIELDINTS_H_
#define FIELDINTS_H_

#include "fieldint.h"

class fieldint1 : public fieldint
{
public:
    fieldint1 (lexer* p) : fieldint(p) {};
    virtual ~fieldint1() = default;
};

class fieldint2 : public fieldint
{
public:
    fieldint2 (lexer* p) : fieldint(p) {};
    virtual ~fieldint2() = default;
};

class fieldint3 : public fieldint
{
public:
    fieldint3 (lexer* p) : fieldint(p) {};
    virtual ~fieldint3() = default;
};

class fieldint4 : public fieldint
{
public:
    fieldint4 (lexer* p) : fieldint(p) {};
    virtual ~fieldint4() = default;
};

#endif
#endif
