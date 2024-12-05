/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef REINITOPO_RK3_H_
#define REINITOPO_RK3_H_

#include"reinitopo.h"
#include"field4a.h"
#include"increment.h"

class reinidisc;
class picard;

class reinitopo_RK3 : public reinitopo, public increment
{
public:
	reinitopo_RK3(lexer* p);
	virtual ~reinitopo_RK3() = default;
	virtual void start(lexer*,fdm*,ghostcell*,field&);

	field4a f,frk1,frk2,L,dt;

private:
	reinidisc *prdisc;

	void step(lexer*, fdm*);
    void time_preproc(lexer*);
	
	double starttime,endtime;

	int gcval,gcval_topo,gcval_initopo,reiniter,n;
	const double epsi;
};

#endif
