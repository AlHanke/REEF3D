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

#include"freesurface.h"
#include"gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class picard;
class heat;
class fluid_update;

using namespace std;

#ifndef VOF_IM1_H_
#define VOF_IM1_H_

class VOF_IM1 : public freesurface, gradient
{
public:
	VOF_IM1(lexer*, fdm*, ghostcell*,heat*);
	virtual ~VOF_IM1();
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particlecorr*,field&);
	virtual void ltimesave(lexer*,fdm*,field&);
	virtual void update(lexer*,fdm*,ghostcell*,field&);
	void timesource(lexer*,fdm*,ioflow*);


	void compression(lexer*,fdm*,ghostcell*,convection*,field&,double);

private:
    fluid_update *pupdate;
	field4 phin;
	field1 uc;
	field2 vc;
	field3 wc;
    field4 F;

	int gcval_frac,count,q;
	double starttime;
	void clearrhs(lexer*,fdm*);
	
	convection *ppconvec;
};
#endif

