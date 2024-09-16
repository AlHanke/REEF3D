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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4a.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void partres::debug(lexer *p, fdm &a, ghostcell &pgc, particles_obj2 &PP, sediment_fdm &s)
{

    particles_obj2 Send[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};
    particles_obj2 Recv[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};

    if(p->mpirank==0)
    cout<<"Sending:"<<endl;
    pgc.gcsync();

    int scaler = p->mpirank+1;
    PP.particles[0].xrk1=scaler;
    PP.particles[0].yrk1=scaler+0.1;
    PP.particles[0].zrk1=scaler+0.2;
    PP.particles[0].u=scaler+0.3;
    PP.particles[0].v=scaler+0.4;
    PP.particles[0].w=scaler+0.5;
    PP.particles[0].urk1=scaler+0.6;
    PP.particles[0].vrk1=scaler+0.7;
    PP.particles[0].wrk1=scaler+0.8;

    Send[3].add_entry(PP.particles[0]);
    cout<<p->mpirank+1<<"\n";
    // PP.print(PP.particles[0]);
    Send[3].print(Send[3].particles[0]);
    cout<<flush;

    pgc.gcsync();

    pgc.para_particles_obj2(p,Send,Recv);

    pgc.gcsync();
    if(p->mpirank==0)
    cout<<"Reciving:"<<endl;
    pgc.gcsync();

    cout<<"\n"<<p->mpirank+1<<"\n";
    // for (int n=0;n<6;n++)
    Recv[0].print(Recv[0].particles[0]);
    cout<<flush;
}