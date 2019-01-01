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

#include"rans_io.h"
#include"bc_ikepsilon.h"

class vrans;

using namespace std;

#ifndef IKEPSILON_H_
#define IKEPSILON_H_

class ikepsilon : public rans_io, public bc_ikepsilon
{
public:
	ikepsilon(lexer*,fdm*,ghostcell*);
	virtual ~ikepsilon();
	virtual void isource(lexer*,fdm*);
	virtual void jsource(lexer*,fdm*);
	virtual void ksource(lexer*,fdm*);
	virtual void kinsource(lexer*,fdm*);
	virtual void epssource(lexer*,fdm*);
	virtual void epsfsf(lexer*,fdm*,ghostcell*);
	virtual void eddyvisc(fdm*,lexer*,ghostcell*);
	virtual void clearfield(lexer*,fdm*,field&);

	int count,q;
	double starttime;
    
    vrans *pvrans;
};

#endif



