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
#include"sliceint.h"

void ghostcell::gcsl_neumann_int(sliceint &f, int cs)
{
    for(q=0;q<margin;++q)
    {
        if(cs==1)
            f(i-q-1,j)=f(i,j);
        else if(cs==2)
            f(i,j+q+1)=f(i,j);
        else if(cs==3)
            f(i,j-q-1)=f(i,j);
        else if(cs==4)
            f(i+q+1,j)=f(i,j);
    }
}

void ghostcell::gcsl_neumann_V_int(int *f, int cs)
{
    for(q=0;q<margin;++q)
    {
        if(cs==1)
        {
            f[Im1J]=f[IJ];
            f[Im2J]=f[IJ];
            f[Im3J]=f[IJ];
        }
        else if(cs==2)
        {
            f[IJp1]=f[IJ];
            f[IJp2]=f[IJ];
            f[IJp3]=f[IJ];
        }
        else if(cs==3)
        {
            f[IJm1]=f[IJ];
            f[IJm2]=f[IJ];
            f[IJm3]=f[IJ];
        }
        else if(cs==4)
        {
            f[Ip1J]=f[IJ];
            f[Ip2J]=f[IJ];
            f[Ip3J]=f[IJ];
        }
    }
}
