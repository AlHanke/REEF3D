/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"multiphase_v.h"

void multiphase_v::ini(lexer*, fdm*, ghostcell*, ioflow*, printer*, convection*, solver*)
{
}

void multiphase_v::start(lexer*, fdm*, ghostcell*, convection*, solver*, ioflow*, reini*, particle_corr*, printer*)
{
}

void multiphase_v::update(lexer*, fdm*, ghostcell*)
{
}

void multiphase_v::print_3D(lexer*, fdm*, ghostcell*,ofstream&)
{
}

void multiphase_v::print_file(lexer*, fdm*, ghostcell*)
{
}

double multiphase_v::ls1val(int,int,int)
{
    return 0.0;
}

double multiphase_v::ls2val(int,int,int)
{
    return 0.0;
}

double multiphase_v::ccipol_ls1val(lexer*,ghostcell*,double,double,double)
{
    return 0.0;
}

double multiphase_v::ccipol_ls2val(lexer*,ghostcell*,double,double,double)
{
    return 0.0;
}

void multiphase_v::ls1get(int,int,int,double)
{
}

void multiphase_v::ls2get(int,int,int,double)
{
}

void multiphase_v::name_ParaView_parallel(lexer*, fdm*, ghostcell*,ofstream&)
{
}

void multiphase_v::name_ParaView(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)
{
}

void multiphase_v::offset_ParaView(lexer*, int*, int &)
{
}
