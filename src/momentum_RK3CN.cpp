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
Authors: Elyas Larkermani, Hans Bihs
--------------------------------------------------------------------*/

#include"momentum_RK3CN.h"
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
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"

momentum_RK3CN::momentum_RK3CN(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver,
                                                    ioflow *pioflow, fsi *ppfsi)
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

    if(p->W90>0 && p->F300==0)
    pupdate = new fluid_update_rheology(p);
    else
    pupdate = new fluid_update_void();
}

momentum_RK3CN::~momentum_RK3CN()
{
    delete pupdate;
}

void momentum_RK3CN::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
    pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
    pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);

//Step 1
//--------------------------------------------------------

    // U
    double starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);

    irhs(p,a);
    pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    addirhs(p,a,1.0);
    pdiff->diff_u(p,a,pgc,psolv,urk1,a->u,a->u,a->v,a->w,8.0/15.0);

    p->utime=pgc->timer()-starttime;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);

    jrhs(p,a);
    pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
    addjrhs(p,a,1.0);
    pdiff->diff_v(p,a,pgc,psolv,vrk1,a->v,a->u,a->v,a->w,8.0/15.0);

    p->vtime=pgc->timer()-starttime;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);

    krhs(p,a);
    pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
    addkrhs(p,a,1.0);
    pdiff->diff_w(p,a,pgc,psolv,wrk1,a->w,a->u,a->v,a->w,8.0/15.0);

    p->wtime=pgc->timer()-starttime;

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 8.0/15.0);

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

    pupdate->start(p,a,pgc,a->u,a->v,a->w);

//Step 2
//--------------------------------------------------------

    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);

    irhs(p,a);
    addirhs(p,a,1.0);
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    addirhs(p,a,25.0/8.0);
    pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    addirhs(p,a,-17.0/8.0);
    pdiff->diff_u(p,a,pgc,psolv,urk2,urk1,urk1,vrk1,wrk1,2.0/15.0);

    p->utime+=pgc->timer()-starttime;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);

    jrhs(p,a);
    addjrhs(p,a,1.0);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    addjrhs(p,a,25.0/8.0);
    pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
    addjrhs(p,a,-17.0/8.0);
    pdiff->diff_v(p,a,pgc,psolv,vrk2,vrk1,urk1,vrk1,wrk1,2.0/15.0);

    p->vtime+=pgc->timer()-starttime;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);

    krhs(p,a);
    addkrhs(p,a,1.0);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
    addkrhs(p,a,25.0/8.0);
    pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
    addkrhs(p,a,-17/8.0);
    pdiff->diff_w(p,a,pgc,psolv,wrk2,wrk1,urk1,vrk1,wrk1,2.0/15.0);

    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 2.0/15.0);

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

    pupdate->start(p,a,pgc,a->u,a->v,a->w);

//Step 3
//--------------------------------------------------------

    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);

    irhs(p,a);
    addirhs(p,a,1.0);
    pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    addirhs(p,a,9.0/4.0);
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    addirhs(p,a,-5.0/4.0);
    pdiff->diff_u(p,a,pgc,psolv,a->u,urk2,urk2,vrk2,wrk2,1.0/3.0);

    p->utime+=pgc->timer()-starttime;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);

    jrhs(p,a);
    addjrhs(p,a,1.0);
    pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
    addjrhs(p,a,9.0/4.0);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    addjrhs(p,a,-5.0/4.0);
    pdiff->diff_v(p,a,pgc,psolv,a->v,vrk2,urk2,vrk2,wrk2,1.0/3.0);

    p->vtime+=pgc->timer()-starttime;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);

    krhs(p,a);
    addkrhs(p,a,1.0);
    pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
    addkrhs(p,a,9.0/4.0);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
    addkrhs(p,a,-5.0/4.0);
    pdiff->diff_w(p,a,pgc,psolv,a->w,wrk2,urk2,vrk2,wrk2,1.0/3.0);

    p->wtime+=pgc->timer()-starttime;

    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,1.0/3.0);

    pflow->u_relax(p,a,pgc,a->u);
    pflow->v_relax(p,a,pgc,a->v);
    pflow->w_relax(p,a,pgc,a->w);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    #else
    pgc->start1(p, a->u, gcval_u);
    pgc->start2(p, a->v, gcval_v);
    pgc->start3(p, a->w, gcval_w);
    #endif

    pupdate->start(p,a,pgc,a->u,a->v,a->w);
}

inline void momentum_RK3CN::irhs(lexer *p, fdm *a)
{
    double dens;
    n = 0;
    ULOOP
    {
        a->maxF=std::max(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
        if (p->H10>0 && p->W90==0 && p->H3==2)
            dens = ((a->dro(i+1,j,k)+a->dro(i,j,k))/(a->ro(i+1,j,k)+a->ro(i,j,k)));
        else
            dens = 1.0;

        a->F(i,j,k) += (a->rhsvec.V[n] + a->gi*dens + p->W29_x + a->Fext(i,j,k))*PORVAL1;

        a->rhsvec.V[n] = 0.0;
        a->Fext(i,j,k) = 0.0;
        ++n;
    }
}

inline void momentum_RK3CN::addirhs(lexer *p, fdm *a, double alpha)
{
    n = 0;
    ULOOP
    {
        a->rhsvec.V[n] += alpha*a->F(i,j,k);
        a->F(i,j,k) = 0.0;
        ++n;
    }
}

inline void momentum_RK3CN::jrhs(lexer *p, fdm *a)
{
    double dens;
    n = 0;
    VLOOP
    {
        a->maxG=std::max(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
        if (p->H10>0 && p->W90==0 && p->H3==2)
            dens = ((a->dro(i,j+1,k)+a->dro(i,j,k))/(a->ro(i,j+1,k)+a->ro(i,j,k)));
        else
            dens = 1.0;

        a->G(i,j,k) += (a->rhsvec.V[n] + a->gj*dens + p->W29_y + a->Gext(i,j,k))*PORVAL2;

        a->rhsvec.V[n] = 0.0;
        a->Gext(i,j,k) = 0.0;
        ++n;
    }
}

inline void momentum_RK3CN::addjrhs(lexer *p, fdm *a, double alpha)
{
    n = 0;
    VLOOP
    {
       a->rhsvec.V[n] += alpha*a->G(i,j,k);
       a->G(i,j,k) = 0.0;
       ++n;
    }
}

inline void momentum_RK3CN::krhs(lexer *p, fdm *a)
{
    double dens;
    n = 0;
    WLOOP
    {
        a->maxH=std::max(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
        if (p->H10>0 && p->W90==0 && p->H3==2)
            dens = ((a->dro(i,j,k+1)+a->dro(i,j,k))/(a->ro(i,j,k+1)+a->ro(i,j,k)));
        else
            dens = 1.0;

        a->H(i,j,k) += (a->rhsvec.V[n] + a->gk*dens + p->W29_z + a->Hext(i,j,k))*PORVAL3;

        a->rhsvec.V[n] = 0.0;
        a->Hext(i,j,k) = 0.0;
        ++n;
    }
}

inline void momentum_RK3CN::addkrhs(lexer *p, fdm *a, double alpha)
{
    n = 0;
    WLOOP
    {
       a->rhsvec.V[n] += alpha*a->H(i,j,k);
       a->H(i,j,k) = 0.0;
       ++n;
    }
}
