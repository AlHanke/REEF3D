/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   
void sixdof_sflow::ini(lexer *p, ghostcell *pgc)
{
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_shipwave(p, pgc);
}

void sixdof_sflow::initialize(lexer *, fdm *, ghostcell *, vector<net*>& )
{
}

void sixdof_sflow::initialize(lexer *, fdm_nhf *, ghostcell *, vector<net*>& )
{
}


void sixdof_sflow::start_cfd(lexer* , fdm* , ghostcell* , vrans* , vector<net*>& , int , field &, field &, field &, field &, field &, field &, bool )
{
}

void sixdof_sflow::start_nhflow(lexer* , fdm_nhf* , ghostcell* , vrans* , vector<net*>& , int , 
                                        double *, double *, double *, slice &, slice &, bool )
{
}