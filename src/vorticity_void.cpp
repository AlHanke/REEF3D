/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"vorticity_void.h"

vorticity_void::vorticity_void(lexer *p, fdm *a)
{
}

vorticity_void::~vorticity_void()
{
}

void vorticity_void::print_3D(lexer* p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, int &m)
{
}

void vorticity_void::name_ParaView_parallel(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void vorticity_void::name_ParaView(lexer *p, fdm *a, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
}

void vorticity_void::offset_ParaView(lexer *p, int *offset, int &n)
{
}


