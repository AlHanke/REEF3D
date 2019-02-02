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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"

void iowave::fnpf_precalc_dirichlet(lexer *p, ghostcell *pgc)
{
    cout<<p->mpirank<<" FNPF  001"<<endl;
    double fsfloc;
    
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        

        eta(i,j) = wave_eta(p,pgc,x,y);
        etaval[count] = eta(i,j);
        
        z = eta(i,j);
        Fifsfval[count] = wave_u(p,pgc,xg,yg,z);
        ++count;
        }
        
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        
            FKLOOP
            {
            z=p->ZSN[FIJK]-p->phimean;
            
            Fival[count] = wave_u(p,pgc,xg,yg,z);
            ++count;
            }
        }
        
}


