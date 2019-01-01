/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"sflow.h"
#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm2D;
class ghostcell;
class sflow_timestep;
class sflow_momentum;
class sflow_pressure;
class solver2D;
class ioflow;
class sflow_fsf;
class sflow_vtp;
class sflow_vtp_bed;
class sflow_convection;
class sflow_diffusion;
class sflow_boussinesq;
class sflow_filter;

using namespace std;

#ifndef SFLOW_F_H_
#define SFLOW_F_H_

class sflow_f : public sflow, public increment
{
public:
	sflow_f(lexer*, fdm2D*,ghostcell*);
	virtual ~sflow_f();
	
	virtual void start(lexer*, fdm2D*, ghostcell*);
	
private:
	void logic(lexer*, fdm2D*, ghostcell*);
	void ini(lexer*, fdm2D*, ghostcell*);
	void loop(lexer*, fdm2D*, ghostcell*);
	
	void maxcoor(lexer*, fdm2D*, ghostcell*);
	
	void print_debug(lexer*, fdm2D*, ghostcell*);
	
	sflow_timestep *ptime;
	sflow_momentum *pmom;
	sflow_pressure *ppress;
	solver2D *psolv;
	solver2D *ppoissonsolv;
	ioflow *pflow;
	sflow_fsf *pfsf;
	sflow_vtp *pprint;
	sflow_vtp_bed *pprintbed;
	sflow_convection *pconvec;
	sflow_diffusion *pdiff;
    sflow_boussinesq *pbouss;
    sflow_filter *pfilter;
	
	double starttime, endtime;
};

#endif