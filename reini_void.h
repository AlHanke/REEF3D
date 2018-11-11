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

using namespace std;

#ifndef REINI_VOID_H_
#define REINI_VOID_H_

class reini_void : public reini
{
public:
	reini_void(lexer* p);
	virtual ~reini_void();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);
	virtual void update(fdm*, lexer*,ghostcell*);
};

#endif
