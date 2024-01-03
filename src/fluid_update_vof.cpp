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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fluid_update_vof::fluid_update_vof(lexer *p, fdm* a, ghostcell* pgc) : dx(p->DXM),
												visc_air(p->W4),visc_water(p->W2),ro_air(p->W3),ro_water(p->W1),visc_body(p->X44)
{
    gcval_ro=1;
	gcval_visc=1;
}

fluid_update_vof::~fluid_update_vof()
{
}

void fluid_update_vof::start(lexer *p, fdm* a, ghostcell* pgc)
{
	double H=0.0;

}


int fluid_update_vof::iocheck;
int fluid_update_vof::iter;