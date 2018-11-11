/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"fluid_update.h"
#include"increment.h"

class fdm;
class lexer;
class ghostcell;
class heat;

using namespace std;

#ifndef FLUID_UPDATE_FSF_HEAT_H_
#define FLUID_UPDATE_FSF_HEAT_H_

class fluid_update_fsf_heat : public fluid_update, increment
{
public:
    fluid_update_fsf_heat(lexer*, fdm*, ghostcell*, heat*&);
	virtual ~fluid_update_fsf_heat();

	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void start3(lexer*, fdm*, ghostcell*,field&,field&);

private:
    void material(lexer*, fdm*, ghostcell*);
    double material_ipol(double**,int,double);

    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	const double dx;
    double epsi;
	double visc_air,visc_water,ro_air,ro_water,alpha_air,alpha_water;
	double **water_density;
	double **water_viscosity;
	double **air_density;
	double **air_viscosity;

	int water_density_num;
	int water_viscosity_num;
	int air_density_num;
	int air_viscosity_num;

	heat *pheat;

};

#endif


