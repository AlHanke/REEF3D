/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"reduction_deyemp.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

reduction_deyemp::reduction_deyemp(lexer *p) : bedslope(p)
{
}

reduction_deyemp::~reduction_deyemp()
{
}

void reduction_deyemp::start(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double r=1.0;
    

    SLICELOOP4
    {
	s->alpha(i,j) = fabs(s->alpha(i,j));

	r = 0.954*pow(1.0-s->teta(i,j)/s->phi(i,j), 0.745)*pow(1.0-s->alpha(i,j)/s->phi(i,j),0.372);

    // limiter
	if( 1.0-s->teta(i,j)/s->phi(i,j) < 0.0 || 1.0-s->alpha(i,j)/s->phi(i,j)< 0.0)
    {
        if(p->S84==1)
        {
        r = cos(s->teta(i,j))*(1.0 - tan(s->teta(i,j)/tan(s->phi(i,j))));
        r*= cos(s->alpha(i,j))*(1.0 - pow(tan(s->alpha(i,j)),2.0)/pow(tan(s->phi(i,j)),2.0));
        }
        
        if(p->S84==2)
        r = 0.1/(fabs(s->gamma(i,j)) + 0.0000001)+0.1;
    }


    r = MAX(r,0.01);
    r = MIN(r,1.25);

	if(p->pos_x()<p->S71)
	r=1.0;

	if(p->pos_x()>p->S72)
	r=10.0;
	
    s->reduce(i,j)=r;
    }
}



