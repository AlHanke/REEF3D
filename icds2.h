/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"discrete.h"
#include"increment.h"

class flux;

#ifndef ICDS2_H_
#define ICDS2_H_

using namespace std;

class icds2 : public discrete,  public increment
{

public:

	icds2 (lexer *);
	virtual ~icds2();

	virtual void start(lexer*,fdm*,field&,int,field&,field&,field&);

private:
	double dx,dy,dz;
	double L;

	int count,rocount,countN,coliN;
	int *range;
    
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

	void aij(lexer*, fdm*, field&, int,field&,field&,field&);
    
    flux *pflux;
};

#endif
