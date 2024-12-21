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

#ifndef VOF_AB_H_
#define VOF_AB_H_

#include"freesurface.h"
#include"gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class picard;
class convection;
class fluid_update;

class VOF_AB : public freesurface, gradient
{
public:
	VOF_AB(lexer*, fdm*, ghostcell*);
	virtual ~VOF_AB() = default;
	virtual void start(lexer*,fdm*,ghostcell*, convection*, ioflow*, reini*, particle_corr*,field&);
	virtual void update(lexer*,fdm*,ghostcell*);

	void compression(lexer*,fdm*,ghostcell*,convection*,field&,double);
	
private:
    fluid_update *pupdate;
	
	field1 uc;
	field2 vc;
	field3 wc;
    field4 F;
	field4 lab;

	int gcval_frac;
	double starttime;
	
	convection *ppconvec;
};
#endif

