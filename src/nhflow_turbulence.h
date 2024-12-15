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

class fdm_nhf;
class lexer;
class nhflow_scalar_convection;
class nhflow_diffusion;
class solver;
class ghostcell;
class ioflow;
class vrans;

#include<fstream>
#include<vector>
#include<sstream>

#ifndef NHFLOW_TURBULENCE_H_
#define NHFLOW_TURBULENCE_H_

using namespace std;

class nhflow_turbulence
{

public:
	virtual void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans*)=0;
	virtual void ktimesave(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual void etimesave(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual void isource(lexer*, fdm_nhf*)=0;
	virtual void jsource(lexer*, fdm_nhf*)=0;
	virtual void ksource(lexer*,fdm_nhf*)=0;

	virtual void print_3D(lexer*, fdm_nhf*, ghostcell*, std::vector<char>&, int&)=0;
    virtual void ini(lexer*, fdm_nhf*, ghostcell*)=0;
    virtual double kinval(int,int,int)=0;
    virtual double epsval(int,int,int)=0;
	virtual void gcupdate(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual double ccipol_kinval(lexer*,ghostcell*,double,double,double)=0;
	virtual double ccipol_epsval(lexer*,ghostcell*,double,double,double)=0;
    virtual double ccipol_a_kinval(lexer*,ghostcell*,double,double,double)=0;
	virtual double ccipol_a_epsval(lexer*,ghostcell*,double,double,double)=0;
    virtual void kinget(int,int,int,double)=0;
    virtual void epsget(int,int,int,double)=0;

    virtual void name_ParaView_parallel(lexer*, fdm_nhf*, ghostcell*,ofstream&)=0;
    virtual void name_ParaView(lexer*, fdm_nhf*, ghostcell*, stringstream&, int*, int &)=0;
    virtual void offset_ParaView(lexer*, int*, int &)=0;
	
	double uref;
};

#endif
