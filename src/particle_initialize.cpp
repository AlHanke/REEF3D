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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"

void particle_f::ini(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{
    
    // seed
    seed_ini(p,a,pgc);
    gpartnum=pgc->globalisum(partnum);
    allocate(p,a,pgc);
    seed(p,a,pgc);
    make_stationary(p,a,&PP);
    
    // print
    print_vtp(p,a,pgc);
    printcount++;
    gparticle_active = pgc->globalisum(PP.size);
    if(p->mpirank==0)
	{
        cout<<"Particles: active: "<<gparticle_active<<endl;
    }
} 
