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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef INCREMENT_H_
#define INCREMENT_H_

#include "looping.h"

class increment
{
public:
    increment() = default;
    virtual ~increment() = default;
    inline static int i = 0;
    inline static int j = 0;
    inline static int k = 0;
    inline static int n = 0;
    inline static int pip = 0;
    constexpr static int marge = 5;
    inline static int max_i = 0;
    inline static int max_j = 0;
    inline static int max_k = 0;
    inline static int level = 0;
};
#endif
