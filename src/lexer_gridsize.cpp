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

#include"lexer.h"

void lexer::gridsize()
{
    lexer *p = this;
    // -------
    cellnum=0;

    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz+flast; ++k)
    if(flag4[IJK]>0)
        ++cellnum;
    // --------

    Iarray(range_row4, M10+5);
    Iarray(range_col4, M10+5);
    Iarray(range_row7, M10+5);
    Iarray(range_col7, M10+5);
}

void lexer::vellast()
{
    ulast=vlast=wlast=1;

    // parallel boundaries
    if(nb4>=0)
        ulast=0;

    if(nb2>=0)
        vlast=0;

    if(nb6>=0)
        wlast=0;

    // non-parallel perioddic bounbdaries
    if(periodic1>=1)
        ulast=0;

    if(periodic2>=1)
        vlast=0;

    if(periodic3>=1)
        wlast=0;

    flast=0;

    if(A10==3 || A10==5)
        flast=1;
}
