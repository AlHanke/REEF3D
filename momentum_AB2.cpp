/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"momentum_AB2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

momentum_AB2::momentum_AB2(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow)
                                                    :bcmom(p),uab(p),vab(p),wab(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
}

momentum_AB2::~momentum_AB2()
{
}

void momentum_AB2::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{
		
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);


	//--------------------------------------------------------
	//U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	p->utime=pgc->timer()-starttime;

	//--------------------------------------------------------
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	p->vtime=pgc->timer()-starttime;

	//--------------------------------------------------------
	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	p->wtime=pgc->timer()-starttime;
	
	//--------------------------------------------------------
    //U 2
    starttime=pgc->timer();

	if(p->count==1)
	ULOOP
	uab(i,j,k)=a->F(i,j,k);

	ULOOP
	{
	a->u(i,j,k)+=p->dt*CPOR1*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->F(i,j,k) - (p->dt/p->dt_old)*uab(i,j,k));

	uab(i,j,k)=a->F(i,j,k);
	}
	pgc->start1(p,a->u,gcval_u);
	
	p->utime+=pgc->timer()-starttime;

	//--------------------------------------------------------
    //V 2
    starttime=pgc->timer();

	if(p->count==1)
	VLOOP
	vab(i,j,k)=a->G(i,j,k);

	VLOOP
	{
	a->v(i,j,k)+=p->dt*CPOR2*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->G(i,j,k) - (p->dt/p->dt_old)*vab(i,j,k));
	
	vab(i,j,k)=a->G(i,j,k);
	}
	pgc->start2(p,a->v,gcval_v);
	
	p->vtime+=pgc->timer()-starttime;

	//--------------------------------------------------------
    //W 2
    starttime=pgc->timer();

	if(p->count==1)
	WLOOP
	wab(i,j,k)=a->H(i,j,k);

	WLOOP
	{
	a->w(i,j,k)+=p->dt*CPOR3*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->H(i,j,k) - (p->dt/p->dt_old)*wab(i,j,k));								
	
	wab(i,j,k)=a->H(i,j,k);
	}
	pgc->start3(p,a->w,gcval_w);
	
	p->wtime+=pgc->timer()-starttime;

	pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	//--------------------------------------------------------
	// pressure	
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,1.0);

	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	pflow->periodic(a->u,p);
	pflow->periodic(a->v,p);
	pflow->periodic(a->w,p);
	
}

void momentum_AB2::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
    pgc->forcing1(p,a,f,uvel,vvel,wvel,alpha); 
    
	n=0;
	if(p->D20<3)
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n]),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	ULOOP
	{
	a->rhsvec.V[n]+=a->gi;
	++n;
	}
}

void momentum_AB2::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
    pgc->forcing2(p,a,f,uvel,vvel,wvel,alpha);
    
	n=0;
	if(p->D20<3)
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n]),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;	
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	VLOOP
	{
	a->rhsvec.V[n]+=a->gj;
	++n;
	}
}

void momentum_AB2::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
    pgc->forcing3(p,a,f,uvel,vvel,wvel,alpha);
    
	n=0;
	if(p->D20<3)
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n]),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
	a->rhsvec.V[n]=0.0;
	++n;
	}
	
	n=0;
	if(p->D20==3)
	WLOOP
	{
	a->rhsvec.V[n]+=a->gk;
	++n;
	}
}

void momentum_AB2::utimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}

void momentum_AB2::vtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}

void momentum_AB2::wtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}

void momentum_AB2::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_AB2::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_AB2::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}


