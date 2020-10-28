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

#include"momentum_fsi.h"
#include"vrans.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"density_f.h"
#include"ediff2.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"6DOF_fsi.h"
#include"net.h"

momentum_fsi::momentum_fsi
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc, 
    convection *pconvection, 
    diffusion *pdiffusion, 
    pressure* ppressure, 
    poisson* ppoisson,
    turbulence *pturbulence, 
    solver *psolver, 
    solver *ppoissonsolver, 
    ioflow *pioflow
):bcmom(p),urk1(p),vrk1(p),wrk1(p),urk2(p),vrk2(p),wrk2(p),un(p),vn(p),wn(p),uf(p),vf(p),wf(p),gradPx(p),gradPy(p),gradPz(p),flagx(p),flagy(p),flagz(p), flagp(p),fx(p),fy(p),fz(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
	gcval_urk=20;
	gcval_vrk=21;
	gcval_wrk=22;

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    
    
    Xfb = Yfb = Zfb = Kfb = Mfb = Nfb = cd = cq = cl = 0.0;
   
    // Define explicit diffusion for predictor step
	pdiff_e = new ediff2(p);

    pdensity = new density_f(p);
}

momentum_fsi::~momentum_fsi(){}


void momentum_fsi::ini(lexer *p, fdm* a, ghostcell* pgc, sixdof_fsi* p6dof_fsi,vrans* pvrans, vector<net*>& pnet)
{ 
    // Calculate initial forcing term
    ULOOP
    {
        fx(i,j,k) = 0.0; 
    }
        
    VLOOP
    {
       fy(i,j,k) = 0.0;
    }   
    
    WLOOP
    {
        fz(i,j,k) = 0.0;
    }

    //forcing(p, a, pgc, p6dof_fsi, a->u,a->v,a->w,a->u,a->v,a->w,1.0,pvrans,pnet);
}

void momentum_fsi::predictor(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom, vrans *pvrans)
{
}

void momentum_fsi::start(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans)
{	}

void momentum_fsi::starti(lexer* p, fdm* a, ghostcell* pgc, sixdof_fsi* p6dof_fsi, vrans* pvrans, vector<net*>& pnet)
{	
    ULOOP
    {
        un(i,j,k) = a->u(i,j,k);
    }

    VLOOP
    {
        vn(i,j,k) = a->v(i,j,k);
    }

    WLOOP
    {
        wn(i,j,k) = a->w(i,j,k);
    }

    pgc->start1(p,un,gcval_u);
    pgc->start2(p,vn,gcval_v);
    pgc->start3(p,wn,gcval_w);


    // Set inflow 
    double udisctime=0.0;
    double udiscstart=0.0;
    
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
		

//Step 1
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans); 
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,a->u,a->u,a->v,a->w,1.0);
    udiscstart=pgc->timer();
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    udisctime=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	ULOOP
	urk1(i,j,k) = a->u(i,j,k)
				+ p->dt*CPOR1*(a->F(i,j,k));

    p->utime=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,a->v,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	VLOOP
	vrk1(i,j,k) = a->v(i,j,k)
				+ p->dt*CPOR2*(a->G(i,j,k));

    p->vtime=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,a->w,a->u,a->v,a->w,1.0);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);

	WLOOP
	wrk1(i,j,k) = a->w(i,j,k)
				+ p->dt*CPOR3*(a->H(i,j,k));
	
    p->wtime=pgc->timer()-starttime;
   

    forcing(p, a, pgc, p6dof_fsi,urk1,vrk1,wrk1,a->u,a->v,a->w,1.0,pvrans,pnet);
	ULOOP
	urk1(i,j,k) += 1.0*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	vrk1(i,j,k) += 1.0*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	wrk1(i,j,k) += 1.0*p->dt*CPOR3*(fz(i,j,k));
    

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
	
	//urk1.ggcpol(p);
	//vrk1.ggcpol(p);
	//wrk1.ggcpol(p);
	
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, urk1, vrk1, wrk1, 1.0);
	
	pflow->u_relax(p,a,pgc,urk1);
	pflow->v_relax(p,a,pgc,vrk1);
	pflow->w_relax(p,a,pgc,wrk1);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,urk1,gcval_urk);
	pgc->start2(p,vrk1,gcval_vrk);
	pgc->start3(p,wrk1,gcval_wrk);
	
//Step 2
//--------------------------------------------------------
	
	
	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk1,urk1,vrk1,wrk1,0.25);
    udiscstart=pgc->timer();
	pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    udisctime+=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	ULOOP
	urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
				+ 0.25*p->dt*CPOR1*(a->F(i,j,k));
                
    p->utime+=pgc->timer()-starttime;
	
	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
	pdiff->diff_v(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	VLOOP
	vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
				+ 0.25*p->dt*CPOR2*(a->G(i,j,k));
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk1,urk1,vrk1,wrk1,0.25);
	pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
	pdiff->diff_w(p,a,pgc,psolv,urk1,vrk1,wrk1,0.25);

	WLOOP
	wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
				+ 0.25*p->dt*CPOR3*(a->H(i,j,k));

    p->wtime+=pgc->timer()-starttime;


    forcing(p, a, pgc, p6dof_fsi,urk2,vrk2,wrk2,a->u,a->v,a->w,0.25,pvrans,pnet);
	ULOOP
	urk2(i,j,k) += 0.25*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	vrk2(i,j,k) += 0.25*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	wrk2(i,j,k) += 0.25*p->dt*CPOR3*(fz(i,j,k));
	
    pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);
	
	//urk2.ggcpol(p);
	//vrk2.ggcpol(p);
	//wrk2.ggcpol(p);

    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, urk2, vrk2, wrk2, 0.25);
	
	pflow->u_relax(p,a,pgc,urk2);
	pflow->v_relax(p,a,pgc,vrk2);
	pflow->w_relax(p,a,pgc,wrk2);
	pflow->p_relax(p,a,pgc,a->press);
	
	pgc->start1(p,urk2,gcval_urk);
	pgc->start2(p,vrk2,gcval_vrk);
	pgc->start3(p,wrk2,gcval_wrk);


//Step 3
//--------------------------------------------------------

	// U
	starttime=pgc->timer();

	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
	ppress->upgrad(p,a);
	irhs(p,a,pgc,urk2,urk2,vrk2,wrk2,2.0/3.0);
    udiscstart=pgc->timer();
	pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    udisctime+=pgc->timer()-udiscstart;
	pdiff->diff_u(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	ULOOP
	a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR1*(a->F(i,j,k));
	
    p->utime+=pgc->timer()-starttime;

	// V
	starttime=pgc->timer();

	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
	ppress->vpgrad(p,a);
	jrhs(p,a,pgc,vrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
	pdiff->diff_v(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	VLOOP
	a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR2*(a->G(i,j,k));
	
    p->vtime+=pgc->timer()-starttime;

	// W
	starttime=pgc->timer();

	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
	ppress->wpgrad(p,a);
	krhs(p,a,pgc,wrk2,urk2,vrk2,wrk2,2.0/3.0);
	pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
	pdiff->diff_w(p,a,pgc,psolv,urk2,vrk2,wrk2,2.0/3.0);

	WLOOP
	a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
				+ (2.0/3.0)*p->dt*CPOR3*(a->H(i,j,k));
	
    p->wtime+=pgc->timer()-starttime;


    forcing(p, a, pgc, p6dof_fsi,a->u,a->v,a->w,un,vn,wn,2.0/3.0,pvrans,pnet);
	ULOOP
	a->u(i,j,k) += 2.0/3.0*p->dt*CPOR1*(fx(i,j,k));
	VLOOP
	a->v(i,j,k) += 2.0/3.0*p->dt*CPOR2*(fy(i,j,k));
	WLOOP
	a->w(i,j,k) += 2.0/3.0*p->dt*CPOR3*(fz(i,j,k));
	

    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);

	//--------------------------------------------------------
	// pressure
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, a->u, a->v,a->w,2.0/3.0);

    pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

    bool conv = true;
    p6dof_fsi->updateFSI(p,a,pgc,conv);
    p6dof_fsi->print_stl(p,a,pgc);
    p6dof_fsi->print_parameter(p, a, pgc);
    if (p->mpirank == 0)
	{
		cout<<"Ue: "<<p->ufbi<<" Ve: "<<p->vfbi<<" We: "<<p->wfbi<<" Pe: "<<p->pfbi<<" Qe: "<<p->qfbi<<" Re: "<<p->rfbi<<endl;
    }
}

void momentum_fsi::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        ULOOP
        {
            a->maxF = MAX(fabs(a->rhsvec.V[n] + a->gi), a->maxF);
            
            a->F(i,j,k) += (a->rhsvec.V[n] + a->gi)*PORVAL1;
            
            a->rhsvec.V[n] = 0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        ULOOP
        {
            a->rhsvec.V[n] += a->gi;
            
            ++n;
        }
    }
}


void momentum_fsi::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        VLOOP
        {
            a->maxG = MAX(fabs(a->rhsvec.V[n] + a->gj), a->maxG);
            
            a->G(i,j,k) += (a->rhsvec.V[n] + a->gj)*PORVAL2;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        VLOOP
        {
            a->rhsvec.V[n] += a->gj;
            
            ++n;
        }
    }
}

void momentum_fsi::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel, double alpha)
{
	n=0;
	if(p->D20 < 3)
    {
        WLOOP
        {
            a->maxH = MAX(fabs(a->rhsvec.V[n] + a->gk), a->maxH);
            
            a->H(i,j,k) += (a->rhsvec.V[n] + a->gk)*PORVAL3;
            
            a->rhsvec.V[n]=0.0;
            
            ++n;
        }
    }
	
	n=0;
	if(p->D20 == 3)
    {
        WLOOP
        {
            a->rhsvec.V[n] += a->gk;
            
            ++n;
        }
    }
}


void momentum_fsi::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    {
        un(i,j,k) = a->u(i,j,k);
    }
    
    pgc->start1(p,un,gcval_u);
}


void momentum_fsi::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
    VLOOP
    {
        vn(i,j,k) = a->v(i,j,k);
    }

    pgc->start2(p,vn,gcval_v);     
}


void momentum_fsi::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
    WLOOP
    {
        wn(i,j,k) = a->w(i,j,k);
    }
        
    pgc->start3(p,wn,gcval_w);            
}

void momentum_fsi::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_fsi::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}
void momentum_fsi::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv){}


