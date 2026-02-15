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

#include"ghostcell.h"
#include"field.h"

void ghostcell::imagepoint(lexer *p, field& f, double &x_ip, double& val_ip, double dx, int cs)
{
    double y0,y1;
    y1 = 0.0;      // x_j-1
    y0 = f(i,j,k); // x_j

    //fill y[]
    if(cs==dir_labels::X_NEG)
        y1=f(i+1,j,k);
    else if(cs==dir_labels::X_POS)
        y1=f(i-1,j,k);
    else if(cs==dir_labels::Y_NEG)
        y1=f(i,j+1,k);
    else if(cs==dir_labels::Y_POS)
        y1=f(i,j-1,k);
    else if(cs==dir_labels::Z_NEG)
        y1=f(i,j,k+1);
    else if(cs==dir_labels::Z_POS)
        y1=f(i,j,k-1);

    x_ip = -(gamma*dx);

    val_ip = (1.0-gamma)*y0 + gamma*y1;
}
