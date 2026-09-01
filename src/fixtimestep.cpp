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

#include "fixtimestep.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "turbulence.h"
#include <algorithm>
#include <iomanip>
#include <cstdlib>

fixtimestep::fixtimestep(lexer*)
{
}

void fixtimestep::start(fdm* a, lexer* p, ghostcell* pgc, turbulence *pturb)
{
    p->dt=p->dt_old=p->N49;
    p->umax=p->vmax=p->wmax=p->viscmax=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
    p->pressmin=1.0e9;

    // maximum velocities
    ULOOP
    p->umax=std::max(p->umax,fabs(a->u(i,j,k)));

    p->umax=pgc->globalmax(p->umax);

    // TEMPORARY (env REEF_UMAX_LOC): locate the max|u| cell so the residual-velocity driver
    // can be attributed (surface / C-F edge / wall / level). Reports the cell whose |u| equals
    // the global umax, with its level, global index, phi, and lateral-wall/surface context.
    if(std::getenv("REEF_UMAX_LOC"))
    {
        #if USE_AMREX
        double loc=0.0; int llev=-1, li[3]={-1,-1,-1}; double lphi=0.0;
        ULOOP
        {
            const double au=fabs(a->u(i,j,k));
            if(au>loc)
            {
                loc=au; llev=p->level;
                li[0]=i+p->amr_tile_lo.x; li[1]=j+p->amr_tile_lo.y; li[2]=k+p->amr_tile_lo.z;
                lphi=0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
            }
        }
        if(loc>=p->umax*(1.0-1e-9) && loc>0.0)
            std::cout<<"  [umaxloc] |u|="<<loc<<" lev="<<llev
                     <<" ("<<li[0]<<","<<li[1]<<","<<li[2]<<")  phi_face="<<lphi<<std::endl;
        #endif
    }


    VLOOP
    p->vmax=std::max(p->vmax,fabs(a->v(i,j,k)));

    p->vmax=pgc->globalmax(p->vmax);


    WLOOP
    p->wmax=std::max(p->wmax,fabs(a->w(i,j,k)));

    p->wmax=pgc->globalmax(p->wmax);


    if(p->mpirank==0 && (p->count%p->P12==0))
    {
        cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
        cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
        cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }

    // maximum viscosity
    LOOP
    p->viscmax=std::max(p->viscmax, a->visc(i,j,k)+a->eddyv(i,j,k));

    p->viscmax=pgc->globalmax(p->viscmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"viscmax: "<<p->viscmax<<endl;

    //----kin
    LOOP
    p->kinmax=std::max(p->kinmax,pturb->kinval(i,j,k));

    p->kinmax=pgc->globalmax(p->kinmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"kinmax: "<<p->kinmax<<endl;

    //---eps
    LOOP
    p->epsmax=std::max(p->epsmax,pturb->epsval(i,j,k));

    p->epsmax=pgc->globalmax(p->epsmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"epsmax: "<<p->epsmax<<endl;


    //---press
    LOOP
    {
        p->pressmax=std::max(p->pressmax,a->press(i,j,k));
        p->pressmin=std::min(p->pressmin,a->press(i,j,k));
    }

    p->pressmax=pgc->globalmax(p->pressmax);
    p->pressmin=pgc->globalmin(p->pressmin);

}

void fixtimestep::ini(fdm*, lexer* p, ghostcell*)
{
    p->dt=p->dt_old=p->N49;
}
