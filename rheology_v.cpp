/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"rheology_v.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h" a
#include"diff_void.h"
#include"ediff2.h"
#include"idiff2.h"
#include"idiff2_FS.h"

rheology_v::rheology_v(lexer *p, fdm *a) 
{

}

rheology_v::~rheology_v()
{
}

double rheology_v::viscosity(lexer *p, fdm *a, ghostcell *pgc)
{
    val=0.0;
    
    
    return val;
}

void rheology_v::u_source(lexer *p, fdm *a)
{
    
}

void rheology_v::v_source(lexer *p, fdm *a)
{
    
}

void rheology_v::w_source(lexer *p, fdm *a)
{
    
}

void rheology_v::filltau(lexer *p, fdm *a, ghostcell *pgc)
{
}