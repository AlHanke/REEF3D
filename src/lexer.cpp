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

#include "lexer.h"
#include "definitions.h"

lexer::lexer() : coordinates(this), interpolation(this), position(this),
                 flag1(this,DataLocation::FACE_X), flag2(this,DataLocation::FACE_Y), flag3(this,DataLocation::FACE_Z),
                 flag4(this,DataLocation::CELL_CENTERED), flag5(this,DataLocation::CELL_CENTERED),
                 flag7(this,DataLocation::NODE_Z),
                 DF(this),
                 #if USE_AMREX
                 DF1(this, &m_df123, 0), DF2(this, &m_df123, 1), DF3(this, &m_df123, 2),
                 #else
                 DF1(this), DF2(this), DF3(this),
                 #endif
                 IO(this), IOSL(this), flagslice1(this), flagslice2(this), flagslice4(this),
                 mpirank(0)
{
    control::ini_default();

    solveriter=0;

    simtime=0.0;
    poissontime=0.0;
    solidread=toporead=0;
    net_count=0;
    mooring_count=0;
}
