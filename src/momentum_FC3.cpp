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

#include"momentum_FC3.h"
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
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_void.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"
#include <iostream>
#include <iomanip>
#if USE_AMREX
#include"amrex_solver.h"
#endif

// REEF_MOM_PROBE: max|u/v/w| of an RK stage field at named points inside the
// step, to localise which sub-stage first grows the velocity

momentum_FC3::momentum_FC3(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, convection *ppfsfdisc, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow,
                                                    heat *&pheat, concentration *&pconc, reini *ppreini,
                                                    fsi *ppfsi)
                                                    :bcmom(p),
                                                    #if USE_AMREX
                                                    m_rk1(make_mf(p,3,&m_rk1)), m_rk2(make_mf(p,3,&m_rk2)), m_f(make_mf(p,3,&m_f)),
                                                    urk1(p,&m_rk1,0), urk2(p,&m_rk2,0), fx(p,&m_f,0),
                                                    vrk1(p,&m_rk1,1), vrk2(p,&m_rk2,1), fy(p,&m_f,1),
                                                    wrk1(p,&m_rk1,2), wrk2(p,&m_rk2,2), fz(p,&m_f,2),
                                                    #else
                                                    urk1(p),urk2(p),vrk1(p),
                                                    vrk2(p),wrk1(p),wrk2(p),
                                                    fx(p),fy(p),fz(p),
                                                    #endif
                                                    udiff(p),vdiff(p),wdiff(p),
                                                    ls(p),frk1(p),frk2(p)
{
    if(p->F50==1)
    gcval_phi=51;
    else if(p->F50==2)
    gcval_phi=52;
    else if(p->F50==3)
    gcval_phi=53;
    else if(p->F50==4)
    gcval_phi=54;

    pconvec=pconvection;
    pfsfdisc=ppfsfdisc;
    pdiff=pdiffusion;
    ppress=ppressure;
    ppois=ppoisson;
    pturb=pturbulence;
    psolv=psolver;
    ppoissonsolv=ppoissonsolver;
    pflow=pioflow;
    preini=ppreini;
    pfsi=ppfsi;

    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
    pupdate = new fluid_update_fsf(p,a,pgc);

    if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
    pupdate = new fluid_update_fsf_comp(p,a,pgc);

    if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==1)
    pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);

    if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==2)
    pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);

    if(p->F30>0 && p->C10>0 && p->W90==0 && p->F300==0)
    pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);

    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
    pupdate = new fluid_update_rheology(p);

    if(p->F300>0 || p->F30==0)
    pupdate = new fluid_update_void();

    if(p->F46==2)
    ppicard = new picard_f(p);
    else if(p->F46==3)
    ppicard = new picard_lsm(p);
    else
    ppicard = new picard_void(p);

    #if USE_AMREX
    amrex_solve = new amrex_solver(p);
    #endif
}

momentum_FC3::~momentum_FC3()
{
    delete pupdate;
    delete ppicard;
    #if USE_AMREX
    delete amrex_solve;
    #endif
}

void momentum_FC3::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    // REEF_MOM_PROBE: max|u/v/w| of the RK stage field at named points inside the
    // step, to localise which sub-stage first grows the velocity (works in both
    // the AMReX and the plain build, so the two can be diffed directly)
    auto fc3probe = [&](const char* tag, field1 &uu, field2 &vv, field3 &ww)
    {
        if(!std::getenv("REEF_MOM_PROBE")) return;
        double mu=0.0, mv=0.0, mw=0.0;
        ULOOP mu = std::max(mu, fabs(uu(i,j,k)));
        VLOOP mv = std::max(mv, fabs(vv(i,j,k)));
        WLOOP mw = std::max(mw, fabs(ww(i,j,k)));
        mu = pgc->globalmax(mu);
        mv = pgc->globalmax(mv);
        mw = pgc->globalmax(mw);
        if(p->mpirank==0)
        std::cout<<"  [momprobe] "<<tag<<"  |u|="<<mu<<" |v|="<<mv<<" |w|="<<mw<<std::endl;
    };

    const double rk3_total_start = pgc->timer();
    double rk3_setup_time = 0.0;
    double rk3_step1_u_time = 0.0;
    double rk3_step1_v_time = 0.0;
    double rk3_step1_w_time = 0.0;
    double rk3_step1_corr_time = 0.0;
    double rk3_step1_u_source_time = 0.0;
    double rk3_step1_u_bcmom_time = 0.0;
    double rk3_step1_u_pgrad_time = 0.0;
    double rk3_step1_u_rhs_time = 0.0;
    double rk3_step1_u_convec_time = 0.0;
    double rk3_step1_u_diff_time = 0.0;
    double rk3_step1_u_update_time = 0.0;
    double rk3_step1_v_source_time = 0.0;
    double rk3_step1_v_bcmom_time = 0.0;
    double rk3_step1_v_pgrad_time = 0.0;
    double rk3_step1_v_rhs_time = 0.0;
    double rk3_step1_v_convec_time = 0.0;
    double rk3_step1_v_diff_time = 0.0;
    double rk3_step1_v_update_time = 0.0;
    double rk3_step1_w_source_time = 0.0;
    double rk3_step1_w_bcmom_time = 0.0;
    double rk3_step1_w_pgrad_time = 0.0;
    double rk3_step1_w_rhs_time = 0.0;
    double rk3_step1_w_convec_time = 0.0;
    double rk3_step1_w_diff_time = 0.0;
    double rk3_step1_w_update_time = 0.0;
    double rk3_step1_corr_forcing_time = 0.0;
    double rk3_step1_corr_pressure_time = 0.0;
    double rk3_step1_corr_relax_time = 0.0;
    double rk3_step1_corr_gc_time = 0.0;
    double rk3_step2_u_time = 0.0;
    double rk3_step2_v_time = 0.0;
    double rk3_step2_w_time = 0.0;
    double rk3_step2_corr_time = 0.0;
    double rk3_step3_u_time = 0.0;
    double rk3_step3_v_time = 0.0;
    double rk3_step3_w_time = 0.0;
    double rk3_step3_corr_time = 0.0;
    double rk3_step1_fsf_time = 0.0;
    double rk3_step1_fsf_disc_time = 0.0;
    double rk3_step1_fsf_reini_time = 0.0;
    double rk3_step1_fsf_picard_time = 0.0;
    double rk3_step1_fsf_update_time = 0.0;
    double rk3_step2_fsf_time = 0.0;
    double rk3_step2_fsf_disc_time = 0.0;
    double rk3_step2_fsf_reini_time = 0.0;
    double rk3_step2_fsf_picard_time = 0.0;
    double rk3_step2_fsf_update_time = 0.0;
    double rk3_step3_fsf_time = 0.0;
    double rk3_step3_fsf_disc_time = 0.0;
    double rk3_step3_fsf_reini_time = 0.0;
    double rk3_step3_fsf_picard_time = 0.0;
    double rk3_step3_fsf_update_time = 0.0;

    double block_start = pgc->timer();
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
    pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
    pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    rk3_setup_time = pgc->timer() - block_start;

//Step 1
//--------------------------------------------------------

    // FSF
    block_start = pgc->timer();
    a->L.setVal(0.0);

    ls.CopyFrom(a->phi);

    pfsfdisc->start(p,a,ls,4,a->u,a->v,a->w);

    LOOP
    frk1(i,j,k) = ls(i,j,k) + p->dt*a->L(i,j,k);

    // FIELDLOOP(frk1,FIELD_CONST(ls)FIELD_CONST_MEMBER(a,L),
    // frk1(i,j,k) = ls(i,j,k) + p->dt*member_L(i,j,k);
    // )

    pflow->phi_relax(p,pgc,frk1);

    pgc->start4(p,frk1,gcval_phi);
    rk3_step1_fsf_disc_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    p->F44=2;
    preini->start(a,p,frk1, pgc, pflow);
    rk3_step1_fsf_reini_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    ppicard->correct_ls(p,a,pgc,frk1);
    rk3_step1_fsf_picard_time = pgc->timer() - block_start;

    block_start = pgc->timer();

    a->phi.CopyFrom(frk1);

    pgc->start4(p,a->phi,gcval_phi);

    // pupdate->start(p,a,pgc,urk1,vrk1,wrk1);

    pupdate->start(p,a,pgc,a->u,a->v,a->w);
    rk3_step1_fsf_update_time = pgc->timer() - block_start;

    rk3_step1_fsf_time = rk3_step1_fsf_disc_time
                       + rk3_step1_fsf_reini_time
                       + rk3_step1_fsf_picard_time
                       + rk3_step1_fsf_update_time;

    // U
    double starttime=pgc->timer();

    block_start = pgc->timer();
    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    rk3_step1_u_source_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    rk3_step1_u_bcmom_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    #if !USE_AMREX
    ppress->upgrad(p,a,a->eta,a->eta_n);
    #endif
    rk3_step1_u_pgrad_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    irhs(p,a);
    rk3_step1_u_rhs_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
    rk3_step1_u_convec_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pdiff->diff_u(p,a,pgc,psolv,udiff,a->u,a->u,a->v,a->w,1.0);
    rk3_step1_u_diff_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    ULOOP
    urk1(i,j,k) = udiff(i,j,k) + p->dt*CPOR1*a->F(i,j,k);
    rk3_step1_u_update_time = pgc->timer() - block_start;

    p->utime=pgc->timer()-starttime;
    rk3_step1_u_time = p->utime;

    // V
    starttime=pgc->timer();

    block_start = pgc->timer();
    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    rk3_step1_v_source_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    rk3_step1_v_bcmom_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    #if !USE_AMREX
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    #endif
    rk3_step1_v_pgrad_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    jrhs(p,a);
    rk3_step1_v_rhs_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
    rk3_step1_v_convec_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pdiff->diff_v(p,a,pgc,psolv,vdiff,a->v,a->u,a->v,a->w,1.0);
    rk3_step1_v_diff_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    VLOOP
    vrk1(i,j,k) = vdiff(i,j,k) + p->dt*CPOR2*a->G(i,j,k);
    rk3_step1_v_update_time = pgc->timer() - block_start;

    p->vtime=pgc->timer()-starttime;
    rk3_step1_v_time = p->vtime;

    // W
    starttime=pgc->timer();

    block_start = pgc->timer();
    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    rk3_step1_w_source_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    rk3_step1_w_bcmom_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    #if !USE_AMREX
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    #endif
    rk3_step1_w_pgrad_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    krhs(p,a);
    rk3_step1_w_rhs_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
    rk3_step1_w_convec_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pdiff->diff_w(p,a,pgc,psolv,wdiff,a->w,a->u,a->v,a->w,1.0);
    rk3_step1_w_diff_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    WLOOP
    wrk1(i,j,k) = wdiff(i,j,k) + p->dt*CPOR3*a->H(i,j,k);
    rk3_step1_w_update_time = pgc->timer() - block_start;

    p->wtime=pgc->timer()-starttime;
    rk3_step1_w_time = p->wtime;

    fc3probe("s1 diff-only     ",udiff,vdiff,wdiff);
    fc3probe("s1 rhs F/G/H     ",a->F,a->G,a->H);
    fc3probe("s1 post-predictor",urk1,vrk1,wrk1);

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);
    rk3_step1_corr_forcing_time = pgc->timer() - block_start;

    fc3probe("s1 post-forcing  ",urk1,vrk1,wrk1);

    block_start = pgc->timer();
    pflow->pressure_io(p,a,pgc);
    #if USE_AMREX
    amrex_solve->start(p,a,pgc,urk1,vrk1,wrk1,frk1,1.0);
    #else
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
    #endif
    rk3_step1_corr_pressure_time = pgc->timer() - block_start;

    fc3probe("s1 post-project  ",urk1,vrk1,wrk1);

    block_start = pgc->timer();
    pflow->u_relax(p,a,pgc,urk1);
    pflow->v_relax(p,a,pgc,vrk1);
    pflow->w_relax(p,a,pgc,wrk1);
    pflow->p_relax(p,a,pgc,a->press);
    rk3_step1_corr_relax_time = pgc->timer() - block_start;

    fc3probe("s1 post-relax    ",urk1,vrk1,wrk1);

    block_start = pgc->timer();
    #if USE_AMREX
    pgc->startBatch(p, m_rk1, 0, {{&urk1,gcval_u},{&vrk1,gcval_v},{&wrk1,gcval_w}});
    #else
    pgc->start1(p,urk1,gcval_u);
    pgc->start2(p,vrk1,gcval_v);
    pgc->start3(p,wrk1,gcval_w);
    #endif
    rk3_step1_corr_gc_time = pgc->timer() - block_start;

    rk3_step1_corr_time = rk3_step1_corr_forcing_time
                        + rk3_step1_corr_pressure_time
                        + rk3_step1_corr_relax_time
                        + rk3_step1_corr_gc_time;

    pupdate->start(p,a,pgc,urk1,vrk1,wrk1);

//Step 2
//--------------------------------------------------------

    // FSF
    block_start = pgc->timer();
    a->L.setVal(0.0);

    pfsfdisc->start(p,a,frk1,4,urk1,vrk1,wrk1);

    LOOP
    frk2(i,j,k) = 0.75*ls(i,j,k) + 0.25*frk1(i,j,k) + 0.25*p->dt*a->L(i,j,k);

    pflow->phi_relax(p,pgc,frk2);

    pgc->start4(p,frk2,gcval_phi);
    rk3_step2_fsf_disc_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    p->F44=2;
    preini->start(a,p,frk2, pgc, pflow);
    rk3_step2_fsf_reini_time += pgc->timer() - block_start;

    block_start = pgc->timer();
    ppicard->correct_ls(p,a,pgc,frk2);
    rk3_step2_fsf_picard_time += pgc->timer() - block_start;

    block_start = pgc->timer();

    a->phi.CopyFrom(frk2);

    pgc->start4(p,a->phi,gcval_phi);
    rk3_step2_fsf_update_time += pgc->timer() - block_start;

    block_start = pgc->timer();
    // p->F44=2;
    // preini->start(a,p,a->phi, pgc, pflow);
    rk3_step2_fsf_reini_time += pgc->timer() - block_start;

    block_start = pgc->timer();
    // ppicard->correct_ls(p,a,pgc,frk2);
    rk3_step2_fsf_picard_time += pgc->timer() - block_start;

    block_start = pgc->timer();
    pupdate->start(p,a,pgc,urk1,vrk1,wrk1);
    rk3_step2_fsf_update_time += pgc->timer() - block_start;

    rk3_step2_fsf_time = rk3_step2_fsf_disc_time
                       + rk3_step2_fsf_reini_time
                       + rk3_step2_fsf_picard_time
                       + rk3_step2_fsf_update_time;

    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    #if !USE_AMREX
    ppress->upgrad(p,a,a->eta,a->eta_n);
    #endif
    irhs(p,a);
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,urk1,vrk1,wrk1,0.25);

    ULOOP
    urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*udiff(i,j,k) + 0.25*p->dt*CPOR1*a->F(i,j,k);

    p->utime+=pgc->timer()-starttime;
    rk3_step2_u_time = p->utime - rk3_step1_u_time;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    #if !USE_AMREX
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    #endif
    jrhs(p,a);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk1,urk1,vrk1,wrk1,0.25);

    VLOOP
    vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vdiff(i,j,k) + 0.25*p->dt*CPOR2*a->G(i,j,k);

    p->vtime+=pgc->timer()-starttime;
    rk3_step2_v_time = p->vtime - rk3_step1_v_time;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    #if !USE_AMREX
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    #endif
    krhs(p,a);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk1,urk1,vrk1,wrk1,0.25);

    WLOOP
    wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wdiff(i,j,k) + 0.25*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    rk3_step2_w_time = p->wtime - rk3_step1_w_time;

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
    #if USE_AMREX
    amrex_solve->start(p,a,pgc,urk2,vrk2,wrk2,frk2,0.25);
    #else
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);
    #endif

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
    rk3_step2_corr_time = pgc->timer() - block_start;

    pupdate->start(p,a,pgc,urk2,vrk2,wrk2);

//Step 3
//--------------------------------------------------------

    // FSF
    block_start = pgc->timer();
    a->L.setVal(0.0);

    pfsfdisc->start(p,a,frk2,4,urk2,vrk2,wrk2);

    LOOP
    ls(i,j,k) = (1.0/3.0)*ls(i,j,k) + (2.0/3.0)*frk2(i,j,k) + (2.0/3.0)*p->dt*a->L(i,j,k);

    pflow->phi_relax(p,pgc,ls);

    LOOP
    a->phi(i,j,k) =  ls(i,j,k);

    pgc->start4(p,a->phi,gcval_phi);
    rk3_step3_fsf_disc_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    p->F44=3;
    preini->start(a,p,a->phi, pgc, pflow);
    rk3_step3_fsf_reini_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    ppicard->correct_ls(p,a,pgc,a->phi);
    rk3_step3_fsf_picard_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pgc->start4(p,a->phi,gcval_phi);

    pupdate->start(p,a,pgc,urk2,vrk2,wrk2);
    rk3_step3_fsf_update_time = pgc->timer() - block_start;

    rk3_step3_fsf_time = rk3_step3_fsf_disc_time
                       + rk3_step3_fsf_reini_time
                       + rk3_step3_fsf_picard_time
                       + rk3_step3_fsf_update_time;


    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    #if !USE_AMREX
    ppress->upgrad(p,a,a->eta,a->eta_n);
    #endif
    irhs(p,a);
    pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,urk2,vrk2,wrk2,2.0/3.0);

    ULOOP
    a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*udiff(i,j,k) + (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);

    p->utime+=pgc->timer()-starttime;
    rk3_step3_u_time = p->utime - rk3_step1_u_time - rk3_step2_u_time;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    #if !USE_AMREX
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    #endif
    jrhs(p,a);
    pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk2,urk2,vrk2,wrk2,2.0/3.0);

    VLOOP
    a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vdiff(i,j,k) + (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);

    p->vtime+=pgc->timer()-starttime;
    rk3_step3_v_time = p->vtime - rk3_step1_v_time - rk3_step2_v_time;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    #if !USE_AMREX
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    #endif
    krhs(p,a);
    pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk2,urk2,vrk2,wrk2,2.0/3.0);

    WLOOP
    a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wdiff(i,j,k) + (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);

    p->wtime+=pgc->timer()-starttime;
    rk3_step3_w_time = p->wtime - rk3_step1_w_time - rk3_step2_w_time;

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

    pflow->pressure_io(p,a,pgc);
    #if USE_AMREX
    amrex_solve->start(p,a,pgc,a->u,a->v,a->w,a->phi,2.0/3.0);
    #else
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, a->u, a->v, a->w, 2.0/3.0);
    #endif

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
    rk3_step3_corr_time = pgc->timer() - block_start;

    pupdate->start(p,a,pgc,a->u,a->v,a->w);

    const double rk3_total_time = pgc->timer() - rk3_total_start;

    // if(p->mpirank==0)
    // {
    //     const double rk3_measured_time = rk3_setup_time
    //                                + rk3_step1_fsf_time + rk3_step2_fsf_time + rk3_step3_fsf_time
    //                                + rk3_step1_u_time + rk3_step1_v_time + rk3_step1_w_time + rk3_step1_corr_time
    //                                + rk3_step2_u_time + rk3_step2_v_time + rk3_step2_w_time + rk3_step2_corr_time
    //                                + rk3_step3_u_time + rk3_step3_v_time + rk3_step3_w_time + rk3_step3_corr_time;
    //     const double rk3_other_time = rk3_total_time - rk3_measured_time;

    //     const double denom = (rk3_total_time > 0.0) ? rk3_total_time : 1.0;
    //     std::cout<<"momentum FC3 runtime breakdown (s | %total):"<<std::endl;
    //     std::cout<<std::setprecision(6)
    //              <<"  setup: "<<rk3_setup_time<<" | "<<(100.0*rk3_setup_time/denom)<<std::endl
    //              <<"  step1 FSF: "<<rk3_step1_fsf_time<<" | "<<(100.0*rk3_step1_fsf_time/denom)<<std::endl
    //              <<"    step1 FSF disc: "<<rk3_step1_fsf_disc_time<<" | "<<(100.0*rk3_step1_fsf_disc_time/denom)<<std::endl
    //              <<"    step1 FSF reini: "<<rk3_step1_fsf_reini_time<<" | "<<(100.0*rk3_step1_fsf_reini_time/denom)<<std::endl
    //              <<"    step1 FSF picard: "<<rk3_step1_fsf_picard_time<<" | "<<(100.0*rk3_step1_fsf_picard_time/denom)<<std::endl
    //              <<"    step1 FSF update: "<<rk3_step1_fsf_update_time<<" | "<<(100.0*rk3_step1_fsf_update_time/denom)<<std::endl
    //              <<"  step1 U: "<<rk3_step1_u_time<<" | "<<(100.0*rk3_step1_u_time/denom)<<std::endl
    //              <<"  step1 V: "<<rk3_step1_v_time<<" | "<<(100.0*rk3_step1_v_time/denom)<<std::endl
    //              <<"  step1 W: "<<rk3_step1_w_time<<" | "<<(100.0*rk3_step1_w_time/denom)<<std::endl
    //              <<"  step1 corr: "<<rk3_step1_corr_time<<" | "<<(100.0*rk3_step1_corr_time/denom)<<std::endl
    //              <<"    step1 U source: "<<rk3_step1_u_source_time<<" | "<<(100.0*rk3_step1_u_source_time/denom)<<std::endl
    //              <<"    step1 U bcmom: "<<rk3_step1_u_bcmom_time<<" | "<<(100.0*rk3_step1_u_bcmom_time/denom)<<std::endl
    //              <<"    step1 U pgrad: "<<rk3_step1_u_pgrad_time<<" | "<<(100.0*rk3_step1_u_pgrad_time/denom)<<std::endl
    //              <<"    step1 U rhs: "<<rk3_step1_u_rhs_time<<" | "<<(100.0*rk3_step1_u_rhs_time/denom)<<std::endl
    //              <<"    step1 U convec: "<<rk3_step1_u_convec_time<<" | "<<(100.0*rk3_step1_u_convec_time/denom)<<std::endl
    //              <<"    step1 U diff: "<<rk3_step1_u_diff_time<<" | "<<(100.0*rk3_step1_u_diff_time/denom)<<std::endl
    //              <<"    step1 U update: "<<rk3_step1_u_update_time<<" | "<<(100.0*rk3_step1_u_update_time/denom)<<std::endl
    //              <<"    step1 V source: "<<rk3_step1_v_source_time<<" | "<<(100.0*rk3_step1_v_source_time/denom)<<std::endl
    //              <<"    step1 V bcmom: "<<rk3_step1_v_bcmom_time<<" | "<<(100.0*rk3_step1_v_bcmom_time/denom)<<std::endl
    //              <<"    step1 V pgrad: "<<rk3_step1_v_pgrad_time<<" | "<<(100.0*rk3_step1_v_pgrad_time/denom)<<std::endl
    //              <<"    step1 V rhs: "<<rk3_step1_v_rhs_time<<" | "<<(100.0*rk3_step1_v_rhs_time/denom)<<std::endl
    //              <<"    step1 V convec: "<<rk3_step1_v_convec_time<<" | "<<(100.0*rk3_step1_v_convec_time/denom)<<std::endl
    //              <<"    step1 V diff: "<<rk3_step1_v_diff_time<<" | "<<(100.0*rk3_step1_v_diff_time/denom)<<std::endl
    //              <<"    step1 V update: "<<rk3_step1_v_update_time<<" | "<<(100.0*rk3_step1_v_update_time/denom)<<std::endl
    //              <<"    step1 W source: "<<rk3_step1_w_source_time<<" | "<<(100.0*rk3_step1_w_source_time/denom)<<std::endl
    //              <<"    step1 W bcmom: "<<rk3_step1_w_bcmom_time<<" | "<<(100.0*rk3_step1_w_bcmom_time/denom)<<std::endl
    //              <<"    step1 W pgrad: "<<rk3_step1_w_pgrad_time<<" | "<<(100.0*rk3_step1_w_pgrad_time/denom)<<std::endl
    //              <<"    step1 W rhs: "<<rk3_step1_w_rhs_time<<" | "<<(100.0*rk3_step1_w_rhs_time/denom)<<std::endl
    //              <<"    step1 W convec: "<<rk3_step1_w_convec_time<<" | "<<(100.0*rk3_step1_w_convec_time/denom)<<std::endl
    //              <<"    step1 W diff: "<<rk3_step1_w_diff_time<<" | "<<(100.0*rk3_step1_w_diff_time/denom)<<std::endl
    //              <<"    step1 W update: "<<rk3_step1_w_update_time<<" | "<<(100.0*rk3_step1_w_update_time/denom)<<std::endl
    //              <<"    step1 corr forcing: "<<rk3_step1_corr_forcing_time<<" | "<<(100.0*rk3_step1_corr_forcing_time/denom)<<std::endl
    //              <<"    step1 corr pressure: "<<rk3_step1_corr_pressure_time<<" | "<<(100.0*rk3_step1_corr_pressure_time/denom)<<std::endl
    //              <<"    step1 corr relax: "<<rk3_step1_corr_relax_time<<" | "<<(100.0*rk3_step1_corr_relax_time/denom)<<std::endl
    //              <<"    step1 corr gc: "<<rk3_step1_corr_gc_time<<" | "<<(100.0*rk3_step1_corr_gc_time/denom)<<std::endl
    //              <<"  step2 FSF: "<<rk3_step2_fsf_time<<" | "<<(100.0*rk3_step2_fsf_time/denom)<<std::endl
    //              <<"    step2 FSF disc: "<<rk3_step2_fsf_disc_time<<" | "<<(100.0*rk3_step2_fsf_disc_time/denom)<<std::endl
    //              <<"    step2 FSF reini: "<<rk3_step2_fsf_reini_time<<" | "<<(100.0*rk3_step2_fsf_reini_time/denom)<<std::endl
    //              <<"    step2 FSF picard: "<<rk3_step2_fsf_picard_time<<" | "<<(100.0*rk3_step2_fsf_picard_time/denom)<<std::endl
    //              <<"    step2 FSF update: "<<rk3_step2_fsf_update_time<<" | "<<(100.0*rk3_step2_fsf_update_time/denom)<<std::endl
    //              <<"  step2 U: "<<rk3_step2_u_time<<" | "<<(100.0*rk3_step2_u_time/denom)<<std::endl
    //              <<"  step2 V: "<<rk3_step2_v_time<<" | "<<(100.0*rk3_step2_v_time/denom)<<std::endl
    //              <<"  step2 W: "<<rk3_step2_w_time<<" | "<<(100.0*rk3_step2_w_time/denom)<<std::endl
    //              <<"  step2 corr: "<<rk3_step2_corr_time<<" | "<<(100.0*rk3_step2_corr_time/denom)<<std::endl
    //              <<"  step3 FSF: "<<rk3_step3_fsf_time<<" | "<<(100.0*rk3_step3_fsf_time/denom)<<std::endl
    //              <<"    step3 FSF disc: "<<rk3_step3_fsf_disc_time<<" | "<<(100.0*rk3_step3_fsf_disc_time/denom)<<std::endl
    //              <<"    step3 FSF reini: "<<rk3_step3_fsf_reini_time<<" | "<<(100.0*rk3_step3_fsf_reini_time/denom)<<std::endl
    //              <<"    step3 FSF picard: "<<rk3_step3_fsf_picard_time<<" | "<<(100.0*rk3_step3_fsf_picard_time/denom)<<std::endl
    //              <<"    step3 FSF update: "<<rk3_step3_fsf_update_time<<" | "<<(100.0*rk3_step3_fsf_update_time/denom)<<std::endl
    //              <<"  step3 U: "<<rk3_step3_u_time<<" | "<<(100.0*rk3_step3_u_time/denom)<<std::endl
    //              <<"  step3 V: "<<rk3_step3_v_time<<" | "<<(100.0*rk3_step3_v_time/denom)<<std::endl
    //              <<"  step3 W: "<<rk3_step3_w_time<<" | "<<(100.0*rk3_step3_w_time/denom)<<std::endl
    //              <<"  step3 corr: "<<rk3_step3_corr_time<<" | "<<(100.0*rk3_step3_corr_time/denom)<<std::endl
    //              <<"  other: "<<rk3_other_time<<" | "<<(100.0*rk3_other_time/denom)<<std::endl
    //              <<"  total: "<<rk3_total_time<<std::endl;
    // }
}

void momentum_FC3::irhs(lexer *p, fdm *a)
{
    n=0;
    ULOOP
    {
        a->maxF=std::max(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
        a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;

        a->rhsvec.V[n] = 0.0;
        a->Fext(i,j,k) = 0.0;
        ++n;
    }
}

void momentum_FC3::jrhs(lexer *p, fdm *a)
{
    n=0;
    VLOOP
    {
        a->maxG=std::max(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
        a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;

        a->rhsvec.V[n] = 0.0;
        a->Gext(i,j,k) = 0.0;
        ++n;
    }
}

void momentum_FC3::krhs(lexer *p, fdm *a)
{
    n=0;
    WLOOP
    {
        a->maxH=std::max(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
        a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;

        a->rhsvec.V[n] = 0.0;
        a->Hext(i,j,k) = 0.0;
        ++n;
    }
}

