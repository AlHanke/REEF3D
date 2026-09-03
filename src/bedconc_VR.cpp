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

#include "bedconc_VR.h"
#include "lexer.h"
#include "ghostcell.h"
#include "sediment_fdm.h"

#include <algorithm>

bedconc_VR::bedconc_VR(lexer *p) : d50(p->S20)
{
    const double rhosed = p->S22;
    const double rhowat = p->W1;
    const double g = std::sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22);
    const double visc = p->W2;
    const double Rstar = (rhosed-rhowat)/rhowat;

    double Ds = d50*pow((Rstar*g)/(visc*visc),1.0/3.0);

    Ds = Ds>1.0e-10?Ds:1.0e10;

    powDs = pow(Ds,0.3);
}

void bedconc_VR::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double Ts,Tb;
    double Ti,adist;

    // cb* van Rijn
    SEDSLICELOOP
    {
        Ts = s->tau_crit(i,j);
        Tb = s->tau_eff(i,j);

        Ti = std::max((Tb-Ts)/(Ts),0.0);

        if(p->S61==1)
        {
            if(p->A10==5)
            {
                k=0;
                adist = 0.5*p->DZN[KP]*p->WL[IJ];
            }
            else if(p->A10==6)
            {
                k=s->bedk(i,j);
                adist = 0.5*p->DZN[KP];
            }
        }
        else if(p->S61==2)
        {
            adist = 2.0*d50;
        }
        else
        {
            adist = 3.0*d50;
        }

        // powDs = pow(d50*pow((Rstar*g)/(visc*visc),1.0/3.0), 0.3);
        s->cbe(i,j) = std::min((0.015*d50*pow(Ti,1.5))/(powDs*adist), 0.1);
    }

    pgc->gcsl_start4(p,s->cbe,1);
}
