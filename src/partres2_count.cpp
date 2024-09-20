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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres2::count_particles(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    
    particle_count = 0;
    active_count = 0;
    empty_count = 0;
    
    
    
    for(n=0;n<P.index;++n)
    {
    ++particle_count;
    
    if(P.Flag[n]==ACTIVE)
    ++active_count;
    
    if(P.Flag[n]==EMPTY)
    ++empty_count;
    }
    
    
    
    particle_count = pgc->globalsum(particle_count);
    active_count = pgc->globalsum(active_count);
    empty_count = pgc->globalsum(empty_count);
    
    if(p->mpirank==0)
    cout<<"Particle_count: "<<particle_count<<" Active_count: "<<active_count<<" Empty_count: "<<empty_count<<endl;
    
}