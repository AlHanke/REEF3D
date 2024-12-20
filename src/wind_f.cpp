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

#include"wind_f.h"
#include"lexer.h"
#include"slice.h"

wind_f::wind_f(lexer *p) 
{
    wind_forcing_drag_coeff(p);
    
    cosa = cos(p->A571_dir*(PI/180.0));
    sina = sin(p->A571_dir*(PI/180.0));

}

void wind_f::wind_forcing_nhf_x(lexer *p, double *F, slice &WL)
{
    k=p->knoz-1;
    
    SLICELOOP4
    WETDRY
    F[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
}

void wind_f::wind_forcing_nhf_y(lexer *p, double *G, slice &WL)
{
    k=p->knoz-1;
    
    SLICELOOP4
    WETDRY
    G[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*sina;
}

