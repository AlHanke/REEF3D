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

#ifndef SLICE2_H_
#define SLICE2_H_

#include"slice.h"
#include"increment.h"

class slice2 : public slice, increment
{
public:

	slice2 (lexer*);
	virtual ~slice2();

    virtual double& operator()(int, int);
	double& operator[](int);
    virtual void ggcpol(lexer*);
    virtual void resize(lexer*);
    virtual void dealloc(lexer*);
    
	int di,dj;
	int imin,imax,jmax,jmin;

	double ***gcfeld;

private:

	void fieldalloc(lexer *);
	void fieldgcalloc(lexer*);
	void fieldlength(lexer *);

    int iter;
	int gcfeldsize,feldsize;
	
	int rank, gcsl_extra;
	
	double starttime;
	
	lexer *pp;

};

#endif





