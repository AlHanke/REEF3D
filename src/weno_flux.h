/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef WENO_FLUX_H_
#define WENO_FLUX_H_

#include"convection.h"
#include"increment.h"

class flux;

class weno_flux : public convection, public increment
{
public:
    weno_flux(lexer*);
    virtual ~weno_flux() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) final;

private:
    template<typename GenericField>
    inline double aij(lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&, double*, double*, double*);

    template<typename GenericField>
    inline double fx(lexer*, fdm*, const GenericField&, const GenericField&, int, double, int di=0);
    template<typename GenericField>
    inline double fy(lexer*, fdm*, const GenericField&, const GenericField&, int, double, int dj=0);
    template<typename GenericField>
    inline double fz(lexer*, fdm*, const GenericField&, const GenericField&, int, double, int dk=0);

    static constexpr double tttw = 13.0/12.0, fourth=0.4, third=1.0/3.0, sevsix=7.0/6.0, elvsix=11.0/6.0, sixth=1.0/6.0, fivsix=5.0/6.0, tenth=0.1;
    static constexpr double sixten = 0.6, treten = 0.3;
    static constexpr double epsilon = 0.000001;

    flux *pflux;
};

#endif
