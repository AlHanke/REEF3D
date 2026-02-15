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

void ghostcell::gcb_debug(field &f, int cs)
{
    for(q=0; q<margin; ++q)
    {
        if(cs==1)
            f(i-q-1,j,k)=1.0;
        else if(cs==2)
            f(i,j+q+1,k)=1.0;
        else if(cs==3)
            f(i,j-q-1,k)=1.0;
        else if(cs==4)
            f(i+q+1,j,k)=1.0;
        else if(cs==5)
            f(i,j,k-q-1)=1.0;
        else if(cs==6)
            f(i,j,k+q+1)=1.0;
    }
}
