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

#include"benchmark_vortex.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"

benchmark_vortex::benchmark_vortex(lexer *p, fdm *a, ghostcell* pgc)
{
    if(p->global_xmin!=0.0 || p->global_xmax!=1.0 || p->j_dir || p->global_zmin!=0.0 || p->global_zmax!=1.0)
    {
        if(p->mpirank==0)
        std::cerr<<"Warning: benchmark_vortex2D is designed for a unit square domain [0,1]^2. Current domain is ["<<p->global_xmin<<","<<p->global_xmax<<"] x ["<<p->global_zmin<<","<<p->global_zmax<<"].\nAbborting..."<<std::endl;
        pgc->final(true);
    }

    const double xc = 0.5;
    const double zc = 0.75;
    const double radius = 0.15;

    a->phi.setVal(-1.0);
    a->vof.setVal(0.0);

    LOOP
    {
        const double dx = p->pos_x() - xc;
        const double dz = p->pos_z() - zc;
        const double dist = sqrt(dx*dx + dz*dz);
        double sign = dist>radius ? -1.0 : 1.0; // Inside the radius is positive, outside is negative

        a->phi(i,j,k)=sign*dist;

        if(sign==1.0)
        {
            a->vof(i,j,k) = 1.0;
        }
    }
}

void benchmark_vortex::start(lexer* p, fdm *a, ghostcell *pgc, convection *pconvec )
{
    if(p->simtime >= 6.0)
    {
        if(p->mpirank==0)
        std::cout<<"benchmark_vortex2D has completed."<<std::endl;
        pgc->final(false);
    }

    LOOP
    {
        if (p->simtime < 3.0)
        {
            a->u(i,j,k) = -2.0*cos(PI*p->pos_z())*pow(sin(PI*p->pos_x()),2)*sin(PI*p->pos_z());
            a->v(i,j,k) = 0.0;
            a->w(i,j,k) = 2.0*cos(PI*p->pos_x())*pow(sin(PI*p->pos_z()),2)*sin(PI*p->pos_x());
        }
        else if(p->simtime < 6.0)
        {
            a->u(i,j,k) = 2.0*cos(PI*p->pos_z())*pow(sin(PI*p->pos_x()),2)*sin(PI*p->pos_z());
            a->v(i,j,k) = 0.0;
            a->w(i,j,k) = -2.0*cos(PI*p->pos_x())*pow(sin(PI*p->pos_z()),2)*sin(PI*p->pos_x());
        }
    }

    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
    pgc->start3(p,a->w,12);
}
