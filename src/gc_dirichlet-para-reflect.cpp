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
#include"lexer.h"
#include"field.h"

void ghostcell::dirichlet_para_reflect(field& f, double dist, int cs)
{
    double dx;
    if(cs==dir_labels::X_NEG || cs==dir_labels::X_POS)
        dx = p->DXP[IP];
    else if(cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG)
        dx = p->DYP[JP];
    else if(cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS)
        dx = p->DZP[KP];

    int ys=1;
    if(dist>dx*(1.0-1.0e-9) && dist<dx*(1.0+1.0e-9))
        ys=0;

    //fill y[]
    double y[15] = {0.0};
    for(m=0; m<=orderdir-1; m++)
    {
        if(cs==dir_labels::X_NEG )
            y[m] = f(i+orderdir-m-1,j,k);
        else if(cs==dir_labels::X_POS)
            y[m] = f(i-orderdir+m+1,j,k);
        else if(cs==dir_labels::Y_NEG)
            y[m] = f(i,j+orderdir-m-1,k);
        else if(cs==dir_labels::Y_POS)
            y[m] = f(i,j-orderdir+m+1,k);
        else if(cs==dir_labels::Z_NEG)
            y[m] = f(i,j,k+orderdir-m-1);
        else if(cs==dir_labels::Z_POS)
            y[m] = f(i,j,k-orderdir+m+1);
    }

    y[orderdir]=0.0;

    if(ys==1 && dist<gamma*dx)
    {
        double y1;
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

        y[orderdir-2] = (1.0-gamma)*f(i,j,k) + gamma*y1;
    }

    double weight;
    for(q=0; q<margin; ++q)
    {
        y[orderdir+q+1]=0.0;

        for(m=0; m<orderdir; m++)
        {
            weight=0.0;
            for(n=0;n<orderdir;++n)
                if(m==n && q+m==2)
                    weight = -1.0;
            y[orderdir+q+1]+=weight*y[m];
        }
    }

    // write extrapolated ghost cell values into f()
    for(q=0; q<margin; ++q)
    {
        if(cs==dir_labels::X_NEG)
            f(i-q-1,j,k) = y[orderdir+q+ys];
        else if(cs==dir_labels::X_POS)
            f(i+q+1,j,k) = y[orderdir+q+ys];
        else if(cs==dir_labels::Y_NEG)
            f(i,j-q-1,k) = y[orderdir+q+ys];
        else if(cs==dir_labels::Y_POS)
            f(i,j+q+1,k) = y[orderdir+q+ys];
        else if(cs==dir_labels::Z_NEG)
            f(i,j,k-q-1) = y[orderdir+q+ys];
        else if(cs==dir_labels::Z_POS)
            f(i,j,k+q+1) = y[orderdir+q+ys];
    }
}
