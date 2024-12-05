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

#ifndef BCMOM_H_
#define BCMOM_H_

#include"surftens.h"
#include"roughness.h"
class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;

class bcmom : public surftens, public roughness
{
public:
	bcmom(lexer*);
	virtual ~bcmom() = default;
	virtual void bcmom_start(lexer*,fdm*,ghostcell*,turbulence*,field&, int);
	void wall_law_u(lexer*,fdm*,turbulence*,field&,int,int,int,int,int,double);
	void wall_law_v(lexer*,fdm*,turbulence*,field&,int,int,int,int,int,double);
	void wall_law_w(lexer*,fdm*,turbulence*,field&,int,int,int,int,int,double);

private:
	const double kappa;
	double uplus,ks_plus,dist,ks,ustar;
	int ii,jj,kk;
	double value;
	int gcval_phi, bckin;
};
#endif
