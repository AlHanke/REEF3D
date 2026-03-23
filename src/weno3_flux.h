/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef WENO3_FLUX_H_
#define WENO3_FLUX_H_

#include"convection.h"
#include"weno3_nug_func.h"

class flux;

class weno3_flux final : public convection, public weno3_nug_func
{
public:
    weno3_flux(lexer*);
    virtual ~weno3_flux() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    double aij(lexer*, fdm*, field&, int,field&,field&,field&,double*,double*,double*);

    double fx(lexer*, field&, double);
    double fy(lexer*, field&, double);
    double fz(lexer*, field&, double);

    double L,grad;

    double gradx, grady, gradz;
    double fu1,fv1,fw1,fu2,fv2,fw2;

    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    flux *pflux;
};

#endif
