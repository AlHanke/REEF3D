/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

class fdm;
class lexer;
class ghostcell;
class field;

using namespace std;

#ifndef FLUID_UPDATE_H_
#define FLUID_UPDATE_H_

class fluid_update
{
public:

	virtual void start(lexer*, fdm*, ghostcell*)=0;
	virtual void start3(lexer*, fdm*, ghostcell*,field&,field&)=0;

};

#endif

