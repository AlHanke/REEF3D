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

void lexer::flagini2D()
{
    control_calc();

    sliceflagini();

    grid2Dsize();
}

void lexer::gridini2D()
{
    IOSL.resize();
}

void lexer::sliceflagini()
{
    flagslice1.resize();
    flagslice2.resize();
    flagslice4.resize();

    lexer *p = this;

    p->level = 0;
    // TILE_LOOP
    IMALOOP
    JMALOOP
    {
        flagslice4(i,j)=flagslice4_grid[IJ];
        flagslice1(i,j)=flagslice4(i,j);
        flagslice2(i,j)=flagslice4(i,j);
    }
    flagslice4_grid.reset(); // transported into flagslice4; not needed afterwards

    p->level = 0;
    TILE_LOOP
    ILOOP
    JLOOP
    PSLICECHECK4
    {
        if(flagslice4(i+1,j)<0)
        flagslice1(i,j)=flagslice4(i+1,j);

        if(flagslice4(i,j+1)<0)
        flagslice2(i,j)=flagslice4(i,j+1);
    }

    #if USE_AMREX
    flagslice1.fillHigherLevels();
    flagslice2.fillHigherLevels();
    flagslice4.fillHigherLevels();
    #endif

    wet.resize();
}
