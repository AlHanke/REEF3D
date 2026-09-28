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

#include "density_f.h"
#include "lexer.h"
#include "fdm.h"
#include "heaviside_ls.h"

density_f::density_f(lexer* p)
{ 
}

double density_f::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    // exact mean of H over the segment between the two cell centres (see heaviside_ls_avg)
    const double H = heaviside_ls_avg(a->phi(i,j,k), a->phi(i+aa,j+bb,k+cc), p->psi);

    return p->W1*H + p->W3*(1.0-H);
}

bool density_f::roface_segment(lexer *p, double phi_a, double phi_b, double &rho)
{
    const double H = heaviside_ls_avg(phi_a, phi_b, p->psi);
    rho = p->W1*H + p->W3*(1.0-H);
    return true;
}
