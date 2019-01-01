/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"increment.h"

class fdm;
class lexer;

#ifndef DENSITY_H_
#define DENSITY_H_

using namespace std;

class density : virtual public increment
{

public:
    density(lexer*);
	virtual ~density();

	double roface(lexer*,fdm*,int,int,int);
	double ronode(lexer*,fdm*,int,int,int,int);
	
	double H,roval,phival;
	int ii,jj,kk;
	const double epsi,eps;

};

#endif




