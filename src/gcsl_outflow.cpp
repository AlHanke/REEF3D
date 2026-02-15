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
#include"slice.h"

void ghostcell::gcsl_outflow(slice& f, int cs)
{
    for(q=0;q<margin;++q)
    {
        if(cs==X_NEG)
            f(i-q-1,j)=f(i,j);
        else if(cs==Y_POS)
            f(i,j+q+1)=f(i,j);
        else if(cs==Y_NEG)
            f(i,j-q-1)=f(i,j);
    }

	if(cs==X_POS)
    {
        f(i+1,j) = MAX(0.0, f(i,j));

        f(i+2,j) = f(i+1,j);
        f(i+3,j) = f(i+1,j);
    }
}
