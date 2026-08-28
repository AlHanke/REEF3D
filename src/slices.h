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

#ifndef SLICES_H_
#define SLICES_H_

#include "slice.h"

class slice1 final : public slice
{
public:
    slice1(lexer *p) : slice(p) {};
    virtual ~slice1() = default;
};

class slice2 final : public slice
{
public:
    slice2(lexer *p) : slice(p) {};
    virtual ~slice2() = default;
};

class slice4 final : public slice
{
public:
    slice4(lexer *p) : slice(p) {};
    virtual ~slice4() = default;
};

#endif