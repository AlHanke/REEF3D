/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sediment_f::ini(lexer *p, fdm *a,ghostcell *pgc)
{
	double h;

	ILOOP
    JLOOP
	{
		//h=1.0e20;
		KLOOP
		PBASECHECK
		{
		if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
        h = p->ZP[KP]*fabs(a->topo(i,j,k-1)) + p->ZP[KM1]*fabs(a->topo(i,j,k));
		//h = MIN(h,p->ZN[KP]);	
		}
		
		a->bedzh(i,j)=h;
	}
	
	pgc->gcsl_start4(p,a->bedzh,1);
	
	
	SLICELOOP4
	bedtau(i,j)=0.0;
	
	pgc->gcsl_start4(p,bedtau,1);
}