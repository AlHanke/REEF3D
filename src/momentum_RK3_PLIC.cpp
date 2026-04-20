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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"momentum_RK3_PLIC.h"
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
#include"reini.h"
#include"picard.h"
#include"fluid_update_vof.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"nhflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_vof.h"
#include"VOF_PLIC.h"

momentum_RK3_PLIC::momentum_RK3_PLIC(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc,
                                                    fsi *ppfsi)
                                                    :bcmom(p),
                                                    #if USE_AMREX
                                                    m_rk1(make_mf(p,3)), m_rk2(make_mf(p,3)), m_f(make_mf(p,3)),
                                                    urk1(p,&m_rk1,0), urk2(p,&m_rk2,0), fx(p,&m_f,0),
                                                    vrk1(p,&m_rk1,1), vrk2(p,&m_rk2,1), fy(p,&m_f,1),
                                                    wrk1(p,&m_rk1,2), wrk2(p,&m_rk2,2), fz(p,&m_f,2),
                                                    #else
                                                    urk1(p),urk2(p),vrk1(p),
                                                    vrk2(p),wrk1(p),wrk2(p),
                                                    fx(p),fy(p),fz(p),
                                                    #endif
                                                    udiff(p),vdiff(p),wdiff(p)

{
    pconvec=pconvection;
    pdiff=pdiffusion;
    ppress=ppressure;
    ppois=ppoisson;
    pturb=pturbulence;
    psolv=psolver;
    ppoissonsolv=ppoissonsolver;
    pflow=pioflow;
    pfsi=ppfsi;

    pplic=new VOF_PLIC(p,a,pgc,pheat);
}

momentum_RK3_PLIC::~momentum_RK3_PLIC()
{
    delete pplic;
}

void momentum_RK3_PLIC::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
    pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
    pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);

//********************************************************
//Step 1
//********************************************************

    //-------------------------------------------

    #if USE_AMREX
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    #else
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    #endif

    // advect U
    pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
    pconvec->start(p,a,a->w,3,a->u,a->v,a->w);

    ULOOP
    urk1(i,j,k) = a->u(i,j,k)
                + p->dt*CPOR1*a->F(i,j,k);

    VLOOP
    vrk1(i,j,k) = a->v(i,j,k)
                + p->dt*CPOR2*a->G(i,j,k);

    WLOOP
    wrk1(i,j,k) = a->w(i,j,k)
                + p->dt*CPOR3*a->H(i,j,k);

    // clear_FGH
    clear_FGH(p,a);

    #if USE_AMREX
    pgc->startBatch(p, m_rk1, 0, {{&urk1,gcval_u},{&vrk1,gcval_v},{&wrk1,gcval_w}});
    #else
    pgc->start1(p,urk1,gcval_u);
    pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    #endif

    //-------------------------------------------
    // U
    double starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);
    irhs(p,a);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,a->u,a->v,a->w,1.0);

    ULOOP
    urk1(i,j,k) = udiff(i,j,k)
                + p->dt*CPOR1*a->F(i,j,k);

    p->utime=pgc->timer()-starttime;

    //-------------------------------------------
    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    jrhs(p,a);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk1,a->u,a->v,a->w,1.0);

    VLOOP
    vrk1(i,j,k) = vdiff(i,j,k)
                + p->dt*CPOR2*a->G(i,j,k);

    p->vtime=pgc->timer()-starttime;

    //-------------------------------------------
    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    krhs(p,a);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk1,a->u,a->v,a->w,1.0);

    WLOOP
    wrk1(i,j,k) = wdiff(i,j,k)
                + p->dt*CPOR3*a->H(i,j,k);

    p->wtime=pgc->timer()-starttime;

    #if USE_AMREX
    pgc->startBatch(p, m_rk1, 0, {{&urk1,gcval_u},{&vrk1,gcval_v},{&wrk1,gcval_w}});
    #else
    pgc->start1(p,urk1,gcval_u);
    pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    #endif
    clear_FGH(p,a);

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);

    pflow->u_relax(p,a,pgc,urk1);
    pflow->v_relax(p,a,pgc,vrk1);
    pflow->w_relax(p,a,pgc,wrk1);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    pgc->startBatch(p, m_rk1, 0, {{&urk1,gcval_u},{&vrk1,gcval_v},{&wrk1,gcval_w}});
    #else
    pgc->start1(p,urk1,gcval_u);
    pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    #endif

    clear_FGH(p,a);

//********************************************************
//Step 2
//********************************************************

    // advect U
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);

    ULOOP
    urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*urk1(i,j,k)
                + 0.25*p->dt*CPOR1*a->F(i,j,k);

    VLOOP
    vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vrk1(i,j,k)
                + 0.25*p->dt*CPOR2*a->G(i,j,k);

    WLOOP
    wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wrk1(i,j,k)
                + 0.25*p->dt*CPOR3*a->H(i,j,k);

    #if USE_AMREX
    pgc->startBatch(p, m_rk2, 0, {{&urk2,gcval_u},{&vrk2,gcval_v},{&wrk2,gcval_w}});
    #else
    pgc->start1(p,urk2,gcval_u);
    pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    #endif
    clear_FGH(p,a);

    //-------------------------------------------
    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);
    irhs(p,a);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,urk1,vrk1,wrk1,1.0);

    ULOOP
    urk2(i,j,k) = udiff(i,j,k)
                + 0.25*p->dt*CPOR1*a->F(i,j,k);

    p->utime+=pgc->timer()-starttime;

    //-------------------------------------------
    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    jrhs(p,a);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk2,urk1,vrk1,wrk1,1.0);

    VLOOP
    vrk2(i,j,k) = vdiff(i,j,k)
                + 0.25*p->dt*CPOR2*a->G(i,j,k);

    p->vtime+=pgc->timer()-starttime;

    //-------------------------------------------
    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    krhs(p,a);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk2,urk1,vrk1,wrk1,1.0);

    WLOOP
    wrk2(i,j,k) = wdiff(i,j,k)
                + 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

    #if USE_AMREX
    pgc->startBatch(p, m_rk2, 0, {{&urk2,gcval_u},{&vrk2,gcval_v},{&wrk2,gcval_w}});
    #else
    pgc->start1(p,urk2,gcval_u);
    pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    #endif
    clear_FGH(p,a);

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);

    pflow->u_relax(p,a,pgc,urk2);
    pflow->v_relax(p,a,pgc,vrk2);
    pflow->w_relax(p,a,pgc,wrk2);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    pgc->startBatch(p, m_rk2, 0, {{&urk2,gcval_u},{&vrk2,gcval_v},{&wrk2,gcval_w}});
    #else
    pgc->start1(p,urk2,gcval_u);
    pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    #endif
    clear_FGH(p,a);

//********************************************************
//Step 3
//********************************************************

    //-------------------------------------------
    pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
    pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);

    ULOOP
    a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*urk2(i,j,k)
                + (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);

    VLOOP
    a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vrk2(i,j,k)
                + (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);

    WLOOP
    a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wrk2(i,j,k)
                + (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);

    #if USE_AMREX
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    #else
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    #endif
    clear_FGH(p,a);

    //-------------------------------------------
    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);
    irhs(p,a);
    pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,urk2,vrk2,wrk2,1.0);

    ULOOP
    a->u(i,j,k) = udiff(i,j,k)
                + (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);

    p->utime+=pgc->timer()-starttime;

    //-------------------------------------------
    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    jrhs(p,a);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,a->v,urk2,vrk2,wrk2,1.0);

    VLOOP
    a->v(i,j,k) = vdiff(i,j,k)
                + (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);

    p->vtime+=pgc->timer()-starttime;

    //-------------------------------------------
    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmomPLIC_start(a,p,pgc,pturb,pplic,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    krhs(p,a);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,a->w,urk2,vrk2,wrk2,1.0);

    WLOOP
    a->w(i,j,k) = wdiff(i,j,k)
                + (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;

    #if USE_AMREX
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    #else
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    #endif
    clear_FGH(p,a);

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v,a->w,2.0/3.0);

    pflow->u_relax(p,a,pgc,a->u);
    pflow->v_relax(p,a,pgc,a->v);
    pflow->w_relax(p,a,pgc,a->w);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    #else
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    #endif
    clear_FGH(p,a);
}

void momentum_RK3_PLIC::irhs(lexer *p, fdm *a)
{
    n=0;
    ULOOP
    {
        a->maxF=std::max(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
        a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x)*PORVAL1;
        a->rhsvec.V[n]=0.0;
        ++n;
    }
}

void momentum_RK3_PLIC::jrhs(lexer *p, fdm *a)
{
    n=0;
    VLOOP
    {
        a->maxG=std::max(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
        a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y)*PORVAL2;
        a->rhsvec.V[n]=0.0;
        ++n;
    }
}

void momentum_RK3_PLIC::krhs(lexer *p, fdm *a)
{
    n=0;
    WLOOP
    {
        a->maxH=std::max(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
        a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z)*PORVAL3;
        a->rhsvec.V[n]=0.0;
        ++n;
    }
}

void momentum_RK3_PLIC::clear_FGH(lexer *p, fdm *a)
{
    a->F.setVal(0.0);
    a->G.setVal(0.0);
    a->H.setVal(0.0);
}
