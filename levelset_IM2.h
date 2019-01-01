/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"freesurface.h"
#include"gradient.h"
#include"field4.h"

class picard;
class heat;
class concentration;
class fluid_update;
class ioflow;

using namespace std;

#ifndef LEVELSET_IM2_H_
#define LEVELSET_IM2_H_

class levelset_IM2 : public freesurface, gradient
{
public:
	levelset_IM2(lexer*, fdm*, ghostcell*, heat*&, concentration*&, field&);
	virtual ~levelset_IM2();
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particlecorr*,field&);
	virtual void ltimesave(lexer*,fdm*,field&);
    virtual void update(lexer*,fdm*,ghostcell*,field&);
    
	void timesource(lexer*,fdm*,ioflow*);


	field4 phin,phinn;

private:
    fluid_update *pupdate;
    picard *ppicard;

    void clearrhs(lexer*,fdm*);
	int gcval_phi;
	double starttime;
	double aii;
	int count,q;
};

#endif


