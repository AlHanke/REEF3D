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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::stress_tensor(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double maxTs = 0;
    ALOOP
    {
        Ts(i,j,k) = (1.0/6.0)*PI*pow(P.d50,3.0)*cellSum(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        maxTs = max(maxTs,Ts(i,j,k));
    
    Tau(i,j,k) = Ps*pow(Ts(i,j,k),beta)/std::max(Tc-Ts(i,j,k),epsilon*(1.0-Ts(i,j,k)));
    }
    
    pgc->start4a(p,Tau,1);
    pgc->start4a(p,Ts,1);

    maxTs = pgc->globalmax(maxTs);
    if(p->mpirank==0)
        cout<<"maxTs: "<<maxTs<<endl;
}

