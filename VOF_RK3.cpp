/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"VOF_RK3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"discrete.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

VOF_RK3::VOF_RK3(lexer* p, fdm *a, ghostcell* pgc, heat *pheat):gradient(p),uc(p),vc(p),wc(p),F(p)
{
    if(p->F50==1)
	gcval_frac=71;

	if(p->F50==2)
	gcval_frac=72;

	if(p->F50==3)
	gcval_frac=73;

	if(p->F50==4)
	gcval_frac=74;

	pupdate = new fluid_update_vof(p,a,pgc);
	
	ppdisc = new hric(p);
}

VOF_RK3::~VOF_RK3()
{
}

void VOF_RK3::start(fdm* a,lexer* p, discrete* pdisc,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particlecorr* ppart, field &F)
{
    field4 ark1(p),ark2(p);
    
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkin(p,a,pgc,ark2);
    pflow->fsfrkout(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark2);

// Step 1
    starttime=pgc->timer();

    LOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,a->phi,4,a->u,a->v,a->w);

	LOOP
	ark1(i,j,k) = a->phi(i,j,k)
				+ p->dt*a->L(i,j,k);

    compression(p,a,pgc,pdisc,ark1,1.0);

	pgc->start4(p,ark1,gcval_frac);

// Step 2

    LOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark1,4,a->u,a->v,a->w);

	LOOP
	ark2(i,j,k) = 0.75*a->phi(i,j,k)
				+ 0.25*ark1(i,j,k)
				+ 0.25*p->dt*a->L(i,j,k);

    compression(p,a,pgc,pdisc,ark2,0.25);

	pgc->start4(p,ark2,gcval_frac);

// Step 3

    LOOP
	a->L(i,j,k)=0.0;

	pdisc->start(p,a,ark2,4,a->u,a->v,a->w);


	LOOP
	a->phi(i,j,k) = (1.0/3.0)*a->phi(i,j,k)
				  + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->dt*a->L(i,j,k);

    compression(p,a,pgc,pdisc,a->phi,1.0);
	

	pgc->start4(p,a->phi,gcval_frac);


	pflow->periodic(a->phi,p);

	pupdate->start(p,a,pgc);

	p->lsmtime=pgc->timer()-starttime;
	
	if(p->mpirank==0)
	cout<<"voftime: "<<setprecision(3)<<p->lsmtime<<endl;
}

void VOF_RK3::ltimesave(lexer* p, fdm *a, field &F)
{
}

void VOF_RK3::update(lexer *p, fdm *a, ghostcell *pgc, field &F)
{
    pupdate->start(p,a,pgc);
}

void VOF_RK3::compression(lexer* p, fdm *a, ghostcell *pgc, discrete *pdisc, field &f, double alpha)
{
	double di,dj,dk, dnorm,nx,ny,nz;
    double umax,vmax,wmax;
	double vvel, uvel, wvel, uabs;
    double timestep;
    int iter;

    LOOP
    F(i,j,k)=f(i,j,k)*(1.0-f(i,j,k));

    pgc->start4(p,F,gcval_frac);

// x
	ULOOP
	{
	di = xdx(a,a->phi);
	dj = xdy(a,a->phi);
	dk = xdz(a,a->phi);
	
	dnorm=sqrt(di*di + dj*dj + dk*dk);

    nx=di/(dnorm>1.0e-15?dnorm:1.0e20);
	
	
		pip=1;
        uvel=a->u(i,j,k);
        pip=0;

        pip=2;
        vvel=0.25*(a->v(i,j,k) + a->v(i+1,j,k) + a->v(i,j-1,k) + a->v(i+1,j-1,k));
        pip=0;

        pip=3;
        wvel=0.25*(a->w(i,j,k) + a->w(i+1,j,k) + a->w(i+1,j,k-1) + a->w(i,j,k-1));
        pip=0;
		
	uabs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
	
	uc(i,j,k) = p->F84*uabs * nx;

	}
	
	VLOOP
	{
	di = ydx(a,a->phi);
	dj = ydy(a,a->phi);
	dk = ydz(a,a->phi);
	
	dnorm=sqrt(di*di + dj*dj + dk*dk);

    ny=dj/(dnorm>1.0e-15?dnorm:1.0e20);
	
	
		pip=1;
        uvel=0.25*(a->u(i,j,k) + a->u(i,j+1,k) + a->u(i-1,j,k) + a->u(i-1,j+1,k));
        pip=0;

        pip=2;
        vvel=a->v(i,j,k);
        pip=0;

        pip=3;
        wvel=0.25*(a->w(i,j,k) + a->w(i+1,j,k) + a->w(i+1,j,k-1) + a->w(i,j,k-1));
        pip=0;
		
	uabs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
	
	vc(i,j,k) = p->F84*uabs * ny;

	}
	
	
	WLOOP
	{
	di = zdx(a,a->phi);
	dj = zdy(a,a->phi);
	dk = zdz(a,a->phi);
	
	dnorm=sqrt(di*di + dj*dj + dk*dk);

    nz=dk/(dnorm>1.0e-15?dnorm:1.0e20);
	
	
		pip=1;
        uvel=0.25*(a->u(i,j,k) + a->u(i,j,k+1) + a->u(i-1,j,k) + a->u(i-1,j,k+1));
        pip=0;

        pip=2;
        vvel=0.25*(a->v(i,j,k) + a->v(i,j,k+1) + a->v(i,j-1,k) + a->v(i,j-1,k+1));
        pip=0;

        pip=3;
        wvel=a->w(i,j,k);
        pip=0;
		
	uabs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
	
	wc(i,j,k) = p->F84*uabs * nz;

	}
	
	
    pgc->start1(p,uc,14);
    pgc->start2(p,vc,15);
    pgc->start3(p,wc,16);

    umax=vmax=wmax=0.0;

	ULOOP
    umax=MAX(umax,fabs(uc(i,j,k)));
	
	VLOOP
    vmax=MAX(vmax,fabs(vc(i,j,k)));
	
	WLOOP
    wmax=MAX(wmax,fabs(wc(i,j,k)));
    

    timestep = (0.5*MAX(umax,MAX(vmax,wmax)))/p->dx;

    timestep = pgc->globalmax(timestep);

    if(timestep>=alpha*p->dt)
    iter=1;

    if(timestep<alpha*p->dt)
    {
    iter=int((alpha*p->dt)/timestep);

    timestep = (alpha*p->dt)/(double(iter)+1.0);
    ++iter;
    }

    if(p->mpirank==0)
    cout<<"VOF  dt:"<<timestep<<"  iter: "<<iter<<endl;

    for(int qn=0; qn<iter; ++qn)
    {
    LOOP
	a->L(i,j,k)=0.0;

    ppdisc->start(p,a,F,5,uc,vc,wc);

    LOOP
    f(i,j,k)+=p->dt*a->L(i,j,k);
    }
    

}


