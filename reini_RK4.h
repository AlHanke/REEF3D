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

#include"reini.h"
#include"particlecorr.h"
#include"gradient.h"
#include"vec.h"

class picard;
class reinidisc;

using namespace std;

#ifndef REINI_RK4_H_
#define REINI_RK4_H_

class reini_RK4 : public reini, public increment
{
public:
	reini_RK4(lexer* p, fdm *a);
	virtual ~reini_RK4();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);
	
	int *sizeM;
	vec f,frk1,frk2,frk3,L,dt;

private:
    picard *ppicard;
	reinidisc *prdisc;

	void fsfrkioV(lexer*, fdm*, ghostcell*,vec&);
	void step(lexer*, fdm*);
    void time_preproc(lexer*);
	void inisolid(lexer*, fdm*);
	void disc(lexer *p, fdm*, field&);

	double starttime,endtime;

	int gcval_phi,gcval_ro,gcval_iniphi,reiniter;

};

#endif
