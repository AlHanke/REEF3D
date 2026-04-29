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

#if not USE_AMREX
#ifndef FIELDS_H_
#define FIELDS_H_

#include "field.h"

class field1 : public field
{
public:
    field1(lexer* p) : field(p) {};
    virtual ~field1() = default;
};

class field2 : public field
{
public:
    field2(lexer* p) : field(p) {};
    virtual ~field2() = default;
};

class field3 : public field
{
public:
    field3(lexer* p) : field(p) {};
    virtual ~field3() = default;
};

class field4 : public field
{
public:
    field4(lexer* p) : field(p) {};
    virtual ~field4() = default;
};

#endif
#endif
