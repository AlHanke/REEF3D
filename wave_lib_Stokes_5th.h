/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"wave_lib.h"
#include"wave_lib_parameters.h"
#include"increment.h"

#ifndef WAVE_LIB_STOKES_5TH_H_
#define WAVE_LIB_STOKES_5TH_H_

using namespace std;

class wave_lib_Stokes_5th : public wave_lib, public wave_lib_parameters, public increment
{
public:
    wave_lib_Stokes_5th(lexer*, ghostcell*);
	virtual ~wave_lib_Stokes_5th();
    
    double wave_horzvel(lexer*,double,double,double);
    virtual double wave_horzvel_space_sin(lexer*,double,double,double,int);
    virtual double wave_horzvel_space_cos(lexer*,double,double,double,int);
    virtual double wave_horzvel_time_sin(lexer*,int);
    virtual double wave_horzvel_time_cos(lexer*,int);
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_u_space_sin(lexer*,double,double,double,int);
    virtual double wave_u_space_cos(lexer*,double,double,double,int);
    virtual double wave_u_time_sin(lexer*,int);
    virtual double wave_u_time_cos(lexer*,int);
    
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_v_space_sin(lexer*,double,double,double,int);
    virtual double wave_v_space_cos(lexer*,double,double,double,int);
    virtual double wave_v_time_sin(lexer*,int);
    virtual double wave_v_time_cos(lexer*,int);
    
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_w_space_sin(lexer*,double,double,double,int);
    virtual double wave_w_space_cos(lexer*,double,double,double,int);
    virtual double wave_w_time_sin(lexer*,int);
    virtual double wave_w_time_cos(lexer*,int);
    
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_eta_space_sin(lexer*,double,double,int);
    virtual double wave_eta_space_cos(lexer*,double,double,int);
    virtual double wave_eta_time_sin(lexer*,int);
    virtual double wave_eta_time_cos(lexer*,int);
    
    virtual double wave_fi(lexer*,double,double,double);
    virtual double wave_fi_space_sin(lexer*,double,double,double,int);
    virtual double wave_fi_space_cos(lexer*,double,double,double,int);
    virtual double wave_fi_time_sin(lexer*,int);
    virtual double wave_fi_time_cos(lexer*,int);
    
    
    virtual void parameters(lexer*,ghostcell*);
    
private:
    double a11,a22,a31,a33,a42,a44,a51,a53,a55;
    double b22,b31,b42,b44,b53,b55;
    double e2,e4;
    double singamma,cosgamma;
    double vel,eta,fi,T;
};

#endif
