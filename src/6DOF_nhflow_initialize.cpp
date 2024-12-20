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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_nhflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"vrans.h"
   
void sixdof_nhflow::ini(lexer *, ghostcell *)
{
}

void sixdof_nhflow::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{
    if(p->X10==1 || p->X10==2)
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_nhflow(p, d, pgc, pnet);
    
    if(p->X10==3)
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_shipwave(p, pgc);
}

void sixdof_nhflow::initialize(lexer *, fdm *, ghostcell *, vector<net*>& )
{
}

void sixdof_nhflow::start_sflow(lexer *, ghostcell *, int , slice &, slice &, slice &, slice &, slice &, slice &, slice &, bool )
{
}

void sixdof_nhflow::start_cfd(lexer* , fdm* , ghostcell* , vrans* , vector<net*>& , int , field &, field &, field &, field &, field &, field &, bool )
{
}

