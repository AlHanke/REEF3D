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

#include "ietimestep.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "turbulence.h"

#include <iomanip>

ietimestep::ietimestep(lexer*)
{
}

void ietimestep::start(fdm *a, lexer *p, ghostcell *pgc, turbulence *pturb)
{
    p->umax=p->vmax=p->wmax=p->viscmax=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
    p->pressmin=1.0e9;
    p->dt_old=p->dt;

    p->umax=std::max(p->W11_u,p->umax);
    p->umax=std::max(p->W12_u,p->umax);
    p->umax=std::max(p->W13_u,p->umax);
    p->umax=std::max(p->W14_u,p->umax);
    p->umax=std::max(p->W15_u,p->umax);
    p->umax=std::max(p->W16_u,p->umax);

    p->vmax=std::max(p->W11_v,p->vmax);
    p->vmax=std::max(p->W12_v,p->vmax);
    p->vmax=std::max(p->W13_v,p->vmax);
    p->vmax=std::max(p->W14_v,p->vmax);
    p->vmax=std::max(p->W15_v,p->vmax);
    p->vmax=std::max(p->W16_v,p->vmax);

    p->wmax=std::max(p->W11_w,p->wmax);
    p->wmax=std::max(p->W12_w,p->wmax);
    p->wmax=std::max(p->W13_w,p->wmax);
    p->wmax=std::max(p->W14_w,p->wmax);
    p->wmax=std::max(p->W15_w,p->wmax);
    p->wmax=std::max(p->W16_w,p->wmax);

    // maximum velocities
    ULOOP
    p->umax=std::max(p->umax,fabs(a->u(i,j,k)));

    p->umax=pgc->globalmax(p->umax);


    VLOOP
    p->vmax=std::max(p->vmax,fabs(a->v(i,j,k)));

    p->vmax=pgc->globalmax(p->vmax);


    WLOOP
    p->wmax=std::max(p->wmax,fabs(a->w(i,j,k)));

    p->wmax=pgc->globalmax(p->wmax);


    if(p->mpirank==0 && (p->count%p->P12==0))
    {
        cout<<"umax: "<<setprecision(3)<<p->umax<<" \t utime: "<<p->utime<<endl;
        cout<<"vmax: "<<setprecision(3)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
        cout<<"wmax: "<<setprecision(3)<<p->wmax<<" \t wtime: "<<p->wtime<<endl;
    }

    p->umax=std::max(p->umax,p->ufbmax);
    p->vmax=std::max(p->vmax,p->vfbmax);
    p->wmax=std::max(p->wmax,p->wfbmax);

    // rhs globalmax
    a->maxF=pgc->globalmax(a->maxF);
    a->maxG=pgc->globalmax(a->maxG);
    a->maxH=pgc->globalmax(a->maxH);
    p->fbmax=pgc->globalmax(p->fbmax);
    p->fbmax=pgc->globalmax(p->fbmax);

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


    double cu=1.0e10;
    double cv=1.0e10;
    double cw=1.0e10;
    double cb=1.0e10;

    if(p->N50==1)
    {
        LOOP
        {
            double dx = std::min({p->DXN[IP],p->DYN[JP],p->DZN[KP]});

            cu = std::min(cu, 2.0/((sqrt(p->umax*p->umax + p->vmax*p->vmax + p->wmax*p->wmax))/dx
                                + sqrt((4.0*fabs(std::max({a->maxF,a->maxG,a->maxH})))/dx)));
        }
    }
    else if(p->N50==2)
    {
        LOOP
        {
            cu = std::min(cu, 2.0/((sqrt(p->umax*p->umax))/p->DXN[IP]
                                + sqrt((4.0*fabs(a->maxF))/p->DXN[IP])));

            cv = std::min(cv, 2.0/((sqrt(p->vmax*p->vmax))/p->DYN[JP]
                                + sqrt((4.0*fabs(a->maxG))/p->DYN[JP])));

            cw = std::min(cw, 2.0/((sqrt(p->wmax*p->wmax))/p->DZN[KP]
                                + sqrt((4.0*fabs(a->maxH))/p->DZN[KP])));
        }
    }

    cu = std::min({cu,cv,cw});

    p->dt=p->N47*cu;
    p->dt=pgc->timesync(p->dt);

    // fbdt
    LOOP
    {
        double dx = std::min({p->DXN[IP],p->DYN[JP],p->DZN[KP]});

        cb = std::min(cb, 2.0/sqrt((4.0*fabs(p->fbmax))/dx));
    }

    p->fbdt=p->N47*cb;
    p->fbdt=pgc->timesync(p->fbdt);

    a->maxF=0.0;
    a->maxG=0.0;
    a->maxH=0.0;
}

void ietimestep::ini(fdm* a, lexer* p, ghostcell* pgc)
{
    p->umax=p->vmax=p->wmax=p->viscmax=-1e19;

    p->viscmax = std::max(p->W2,p->W4);

    p->umax=std::max(p->W10,p->umax);

    p->umax=std::max(p->W11_u,p->umax);
    p->umax=std::max(p->W12_u,p->umax);
    p->umax=std::max(p->W13_u,p->umax);
    p->umax=std::max(p->W14_u,p->umax);
    p->umax=std::max(p->W15_u,p->umax);
    p->umax=std::max(p->W16_u,p->umax);

    p->vmax=std::max(p->W11_v,p->vmax);
    p->vmax=std::max(p->W12_v,p->vmax);
    p->vmax=std::max(p->W13_v,p->vmax);
    p->vmax=std::max(p->W14_v,p->vmax);
    p->vmax=std::max(p->W15_v,p->vmax);
    p->vmax=std::max(p->W16_v,p->vmax);

    p->wmax=std::max(p->W11_w,p->wmax);
    p->wmax=std::max(p->W12_w,p->wmax);
    p->wmax=std::max(p->W13_w,p->wmax);
    p->wmax=std::max(p->W14_w,p->wmax);
    p->wmax=std::max(p->W15_w,p->wmax);
    p->wmax=std::max(p->W16_w,p->wmax);

    ULOOP
    p->umax=std::max(p->umax,fabs(a->u(i,j,k)));

    p->umax=pgc->globalmax(p->umax);


    VLOOP
    p->vmax=std::max(p->vmax,fabs(a->v(i,j,k)));

    p->vmax=pgc->globalmax(p->vmax);


    WLOOP
    p->wmax=std::max(p->wmax,fabs(a->w(i,j,k)));

    p->wmax=pgc->globalmax(p->wmax);


    p->umax=std::max(p->umax,2.0*p->ufbmax);
    p->umax=std::max(p->umax,2.0*p->vfbmax);
    p->umax=std::max(p->umax,2.0*p->wfbmax);

    p->umax=std::max(p->umax,2.0*p->X210_u);
    p->umax=std::max(p->umax,2.0*p->X210_v);
    p->umax=std::max(p->umax,2.0*p->X210_w);

    p->umax=std::max(p->umax,2.0);


    double cu=1.0e10;
    double cv=1.0e10;
    double cw=1.0e10;

    if(p->N50==1)
    {
        LOOP
        {
            double dx = std::min({p->DXN[IP],p->DYN[JP],p->DZN[KP]});

            cu = std::min(cu, 2.0/((sqrt(p->umax*p->umax + p->vmax*p->vmax + p->wmax*p->wmax))/dx
                                + sqrt((4.0*fabs(std::max({a->maxF,a->maxG,a->maxH})))/dx)));
        }
    }
    else if(p->N50==2)
    {
        LOOP
        {
            cu = std::min(cu, 2.0/((sqrt(p->umax*p->umax))/p->DXN[IP]
                                + sqrt((4.0*fabs(a->maxF))/p->DXN[IP])));

            cv = std::min(cv, 2.0/((sqrt(p->vmax*p->vmax))/p->DYN[JP]
                                + sqrt((4.0*fabs(a->maxG))/p->DYN[JP])));

            cw = std::min(cw, 2.0/((sqrt(p->wmax*p->wmax))/p->DZN[KP]
                                + sqrt((4.0*fabs(a->maxH))/p->DZN[KP])));
        }
    }

    cu = std::min({cu,cv,cw});

    p->dt=p->N47*cu*0.25;
    p->dt = std::max(p->dt,1.0e-6);

    p->dt=pgc->timesync(p->dt);
    p->dt_old=p->dt;

    a->maxF = fabs(a->gi);
    a->maxG = fabs(a->gj);
    a->maxH = fabs(a->gk);
}
