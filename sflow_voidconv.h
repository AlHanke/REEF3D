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

#include"sflow_discrete.h"
#include"increment.h"

class sflow_flux;
class fdm2D;

using namespace std;

#ifndef SFLOW_VOIDCONV_H_
#define SFLOW_VOIDCONV_H_

class sflow_voidconv : public sflow_discrete, public increment
{
public:
	sflow_voidconv(lexer*);
	virtual ~sflow_voidconv();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);



};

#endif
