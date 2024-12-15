/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef BEDLOAD_DIRECTION_H_
#define BEDLOAD_DIRECTION_H_

#include"increment.h"

class lexer;
class ghostcell;
class sediment_fdm;
class slice;
class sliceint;

class bedload_direction : public increment
{
public:

	virtual void start(lexer*, ghostcell*, sediment_fdm*)=0;


};

#endif


