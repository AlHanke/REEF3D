/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"superbee.h"
#include"lexer.h"
#include"fdm.h"

superbee::superbee (lexer *p)
{
}

superbee::~superbee()
{
}

double superbee::iphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i+q1,j,k)-b(i+q2,j,k));
    r=(b(i+n1,j,k)-b(i+n2,j,k))/(fabs(denom)>1.0e-10?denom:1.0e20);

    phi = MAX(0.0, MAX( MIN(2.0*r,1.0), MIN(r,2.0)));

    return phi;
}

double superbee::jphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i,j+q1,k)-b(i,j+q2,k));
    r=(b(i,j+n1,k)-b(i,j+n2,k))/(fabs(denom)>1.0e-10?denom:1.0e20);

    phi = MAX(0.0, MAX( MIN(2.0*r,1.0), MIN(r,2.0)));

    return phi;
}

double superbee::kphi(field& b,int n1, int n2, int q1, int q2)
{
    denom=(b(i,j,k+q1)-b(i,j,k+q2));
    r=(b(i,j,k+n1)-b(i,j,k+n2))/(fabs(denom)>1.0e-10?denom:1.0e20);

    phi = MAX(0.0, MAX( MIN(2.0*r,1.0), MIN(r,2.0)));

    return phi;
}
