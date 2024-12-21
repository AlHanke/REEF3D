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

#ifndef WAVE_LIB_VOID_H_
#define WAVE_LIB_VOID_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"

class wave_lib_void : public wave_lib_precalc, public wave_lib_parameters, virtual public increment
{
public:
    wave_lib_void(lexer*);
	virtual ~wave_lib_void() = default;

    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_fi(lexer*,double,double,double);

    virtual void parameters(lexer*);
    virtual void wave_prestep(lexer*);
    
private:
    double singamma,cosgamma;
};

#endif
