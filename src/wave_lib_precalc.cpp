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

#include"wave_lib_precalc.h"
#include"lexer.h"

wave_lib_precalc::wave_lib_precalc() 
{ 
    vel=eta=fi=T=0.0;
}

// U -------------------------------------------------------------
double wave_lib_precalc::wave_u_space_sin(lexer *, double, double, double, int)
{
    return vel;
}

double wave_lib_precalc::wave_u_space_cos(lexer *, double, double, double, int)
{
    return vel;
}

double wave_lib_precalc::wave_u_time_sin(lexer *, int)
{
    return T;
}

double wave_lib_precalc::wave_u_time_cos(lexer *, int)
{
    return T;
}


// V -------------------------------------------------------------
double wave_lib_precalc::wave_v_space_sin(lexer *, double, double, double, int)
{
    return vel;
}

double wave_lib_precalc::wave_v_space_cos(lexer *, double, double, double, int)
{ 
    return vel;
}

double wave_lib_precalc::wave_v_time_sin(lexer *, int)
{ 
    return T;
}

double wave_lib_precalc::wave_v_time_cos(lexer *, int)
{
    return T;
}


// W -------------------------------------------------------------
double wave_lib_precalc::wave_w_space_sin(lexer *, double, double, double, int)
{
    return vel;
}

double wave_lib_precalc::wave_w_space_cos(lexer *, double, double, double, int)
{
    return vel;
}

double wave_lib_precalc::wave_w_time_sin(lexer *, int)
{
    return T;
}

double wave_lib_precalc::wave_w_time_cos(lexer *, int)
{
    return T;
}

// ETA -------------------------------------------------------------
double wave_lib_precalc::wave_eta_space_sin(lexer *, double, double, int)
{
    return eta;
}

double wave_lib_precalc::wave_eta_space_cos(lexer *, double, double, int)
{
    return eta;
}

double wave_lib_precalc::wave_eta_time_sin(lexer *, int)
{
    return T;
}

double wave_lib_precalc::wave_eta_time_cos(lexer *, int)
{
    return T;
}

// FI -------------------------------------------------------------
void wave_lib_precalc::wave_fi_precalc_xy_ini(lexer*,int)
{
    
}

void wave_lib_precalc::wave_fi_precalc_xy(lexer*,double,double,int)
{
    
}

void wave_lib_precalc::wave_fi_precalc_n(lexer*)
{
    
}
    
double wave_lib_precalc::wave_fi_space_sin(lexer *, double, double, double, int)
{
    
    return fi;
}

double wave_lib_precalc::wave_fi_space_cos(lexer *, double, double, double, int)
{
    
    return fi;
}

double wave_lib_precalc::wave_fi_time_sin(lexer *, int)
{
    
    return fi;
}

double wave_lib_precalc::wave_fi_time_cos(lexer *, int)
{
    return fi;
}
