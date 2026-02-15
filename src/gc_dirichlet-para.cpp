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

void ghostcell::dirichlet_para(field& f, double dist, int cs)
{
    wallvalue=0.0;
    weight=1.0;

    double dx;
    if(cs==dir_labels::X_NEG || cs==dir_labels::X_POS)
        dx = p->DXP[IP];
    else if(cs==dir_labels::Y_POS || cs==dir_labels::Y_NEG)
        dx = p->DYP[JP];
    else if(cs==dir_labels::Z_NEG || cs==dir_labels::Z_POS)
        dx = p->DZP[KP];

    ys=1;
    if(dist>dx*(1.0-1.0e-9) && dist<dx*(1.0+1.0e-9))
        ys=0;

    //fill pos[]
    for(m=0; m<=orderdir-3; m++)
        pos[m]=-dx*double(orderdir-m-2);

    pos[orderdir-2]=0.0;
    pos[orderdir-1]=dist;

    double x[margin];
    for(m=0; m<margin; m++)
        x[m]=dx*double(m+2-ys);

    //fill y[]
    for(m=0; m<=orderdir-2; m++)
    {
        if(cs==dir_labels::X_NEG )
            y[m] = f(i+orderdir-m-2,j,k);
        else if(cs==dir_labels::X_POS)
            y[m] = f(i-orderdir+m+2,j,k);
        else if(cs==dir_labels::Y_NEG)
            y[m] = f(i,j+orderdir-m-2,k);
        else if(cs==dir_labels::Y_POS)
            y[m] = f(i,j-orderdir+m+2,k);
        else if(cs==dir_labels::Z_NEG)
            y[m] = f(i,j,k+orderdir-m-2);
        else if(cs==dir_labels::Z_POS)
            y[m] = f(i,j,k-orderdir+m+2);
    }

    y[orderdir-1]=wallvalue;

    if(ys==1 && dist<gamma*dx)
    {
        imagepoint(p,f,x_ip,val_ip,dx,cs);

        pos[orderdir-2] = x_ip;
        y[orderdir-2] = val_ip;
    }

    for(q=0; q<margin; ++q)
    {
        y[orderdir+q]=0.0;

        for(m=0; m<orderdir; m++)
        {
            weight=1.0;
            for(n=0;n<orderdir;++n)
            {
                if(m!=n)
                    weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
            }
            y[orderdir+q]+=weight*y[m];
        }
    }

    // write extrapolated ghost cell values into f()
    for(q=0; q<margin; ++q)
    {
        if(cs==dir_labels::X_NEG)
            f(i-q-1,j,k) = y[orderdir+q-1+ys];
        else if(cs==dir_labels::X_POS)
            f(i+q+1,j,k) = y[orderdir+q-1+ys];
        else if(cs==dir_labels::Y_NEG)
            f(i,j-q-1,k) = y[orderdir+q-1+ys];
        else if(cs==dir_labels::Y_POS)
            f(i,j+q+1,k) = y[orderdir+q-1+ys];
        else if(cs==dir_labels::Z_NEG)
            f(i,j,k-q-1) = y[orderdir+q-1+ys];
        else if(cs==dir_labels::Z_POS)
            f(i,j,k+q+1) = y[orderdir+q-1+ys];
    }
}
