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

#ifndef HEAT_PRINT_H_
#define HEAT_PRINT_H_

#include"heat.h"
#include"increment.h"
#include"field4.h"
#include"fluid_update.h"
#include<fstream>

class lexer;
class fdm;
class ghostcell;

using namespace std;

class heat_print : public heat, increment
{
public:
    heat_print(lexer*,fdm*);
	virtual ~heat_print() = default;

    void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, int&);
    void heat_ini(lexer*, fdm*, ghostcell*, heat*);
    virtual double val(int,int,int);

    void name_ParaView_parallel(lexer*, fdm*, ghostcell*,ofstream&);
    void name_ParaView(lexer*, fdm*, ghostcell*, stringstream&, int*, int &);
    void offset_ParaView(lexer*, int*, int &);

    field4 T;
	fluid_update *pupdate;

private:

	float ffn;
	double ddn;
	int n,iin;
};

#endif

