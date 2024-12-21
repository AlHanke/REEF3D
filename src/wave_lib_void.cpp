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

#include"wave_lib_void.h"
#include"lexer.h"

wave_lib_void::wave_lib_void(lexer *p) : wave_lib_parameters(p)
{ 
    parameters(p);
    
    if(p->mpirank==0)
        cout<<"Wave_Lib: no wave specified; "<<endl;
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

// U -------------------------------------------------------------
double wave_lib_void::wave_u(lexer*, double, double, double)
{
    return cosgamma*0.0;
}

// V -------------------------------------------------------------
double wave_lib_void::wave_v(lexer*, double, double, double)
{
    return singamma*0.0;
}

// W -------------------------------------------------------------
double wave_lib_void::wave_w(lexer*, double, double, double)
{
    return 0.0;
}

// ETA -------------------------------------------------------------
double wave_lib_void::wave_eta(lexer *, double , double )
{
    return 0.0;
}

// FI -------------------------------------------------------------
double wave_lib_void::wave_fi(lexer*, double, double, double)
{
    return 0.0;
}

//  -------------------------------------------------------------
void wave_lib_void::parameters(lexer*)
{
}

void wave_lib_void::wave_prestep(lexer*)
{
}
