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

#include"momentum_RK3.h"
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
#include"nhflow.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

#if USE_AMREX
#include <AMReX_MultiFab.H>
// Bug #2 velocity-pattern probe (env REEF_VELPAT): dump w along a z-column through the patch top at
// a fixed (i,j) on the fine level, so the divergence-free CHECKERBOARD in urk2 is directly visible
// (divergence probes are blind to it). Compare flip vs noflip: flip shows w alternating sign cell to
// cell near k=22/23 (patch top HIGH C-F), noflip stays smooth. Strip once bug #2 closes.
static void velpat_probe(lexer* p, fdm* a, field& w, const char* tag)
{
    if(!(std::getenv("REEF_VELPAT") && p->nlevs>1)) return;
    const int lev=1, gi=38, gj=34, k0=13, k1=27;
    auto& w_mf=w.GetMultiFab(lev);
    auto& pr_mf=a->press.GetMultiFab(lev);
    for(amrex::MFIter mfi(w_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& vbx=mfi.validbox();
        // laterally the OWNER of (gi,gj); dump the FULL fab z-range (valid + ghost) so the
        // patch-top C-F ghost band (k>23) is visible -- that is where the checkerboard hides.
        if(gi<vbx.smallEnd(0)||gi>vbx.bigEnd(0)||gj<vbx.smallEnd(1)||gj>vbx.bigEnd(1)) continue;
        const amrex::Box& fbx=mfi.fabbox();
        const auto wa=w_mf.const_array(mfi);
        const auto pa=pr_mf.const_array(mfi);
        const auto Ha=a->H.GetMultiFab(lev).const_array(mfi);   // w-momentum RHS accumulator
        const auto ra=a->ro.GetMultiFab(lev).const_array(mfi);  // cell density (from level-set)
        const auto ga=a->grav_pot.GetMultiFab(lev).const_array(mfi); // W22*pos_z (regrid-rebuilt)
        for(int k=k0;k<=k1;++k)
            if(k>=fbx.smallEnd(2)&&k<=fbx.bigEnd(2))
                std::cout<<"  [velpat] "<<tag<<" rank"<<p->mpirank<<" ("<<gi<<","<<gj<<","<<k<<")"
                         <<(k>vbx.bigEnd(2)||k<vbx.smallEnd(2)?"G":" ")  // G = ghost cell
                         <<" w="<<wa(gi,gj,k)<<" H="<<Ha(gi,gj,k)
                         <<" ro="<<ra(gi,gj,k)<<" gpot="<<ga(gi,gj,k)
                         <<" press="<<pa(gi,gj,k)<<std::endl;
    }
}
#else
static void velpat_probe(lexer*, fdm*, field&, const char*) {}
#endif

momentum_RK3::momentum_RK3(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver,
                                                    ioflow *pioflow, fsi *ppfsi)
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

    if(p->W90==0  || p->F300>0)
    pupdate = new fluid_update_void();

    if(p->W90>0 && p->F300==0)
    pupdate = new fluid_update_rheology(p);
}

momentum_RK3::~momentum_RK3()
{
    delete pupdate;
}

void momentum_RK3::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    // AMR C-F low-face preservation: the post-projection velocity ghost exchange (startBatch
    // below) overwrites the fine LOW C-F normal faces that velcorr corrected, breaking
    // D.u = -L.pcorr at those cells. Save the corrected faces before each exchange and restore
    // them after, so the next stage / free-surface advection sees the consistent velocity.
    // Gated by REEF_CF_LOWFACE_RESTORE during validation; no-op for nlevs<=1 and non-ssamg solvers.
    const bool cf_lowface_restore = (std::getenv("REEF_CF_LOWFACE_RESTORE") != nullptr);

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

    a->test.setVal(0.0,true);

    double block_start = pgc->timer();
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
    pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);
    pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);
    rk3_setup_time = pgc->timer() - block_start;

//Step 1
//--------------------------------------------------------

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
    ppress->upgrad(p,a,a->eta,a->eta_n);
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
    // ULOOP
    // urk1(i,j,k) = udiff(i,j,k)
    //             + p->dt*CPOR1*a->F(i,j,k);
    FIELDLOOP(urk1,
        FIELD_CONST(udiff); FIELD_CONST_MEMBER(a, F); FIELD_CONST_MEMBER(a, porosity),
        const double cpor1 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i+1,j,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        urk1(i,j,k) = udiff(i,j,k) + p->dt*cpor1*member_F(i,j,k);
    )
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
    ppress->vpgrad(p,a,a->eta,a->eta_n);
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
    // VLOOP
    // vrk1(i,j,k) = vdiff(i,j,k)
    //             + p->dt*CPOR2*a->G(i,j,k);
    FIELDLOOP(vrk1,
        FIELD_CONST(vdiff); FIELD_CONST_MEMBER(a, G); FIELD_CONST_MEMBER(a, porosity),
        const double cpor2 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j+1,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        vrk1(i,j,k) = vdiff(i,j,k) + p->dt*cpor2*member_G(i,j,k);
    )
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
    ppress->wpgrad(p,a,a->eta,a->eta_n);
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
    // WLOOP
    // wrk1(i,j,k) = wdiff(i,j,k)
    //             + p->dt*CPOR3*a->H(i,j,k);
    FIELDLOOP(wrk1,
        FIELD_CONST(wdiff); FIELD_CONST_MEMBER(a, H); FIELD_CONST_MEMBER(a, porosity),
        const double cpor3 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j,k+1) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        wrk1(i,j,k) = wdiff(i,j,k) + p->dt*cpor3*member_H(i,j,k);
    )
    rk3_step1_w_update_time = pgc->timer() - block_start;

    p->wtime=pgc->timer()-starttime;
    rk3_step1_w_time = p->wtime;

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk1, vrk1, wrk1, fx, fy, fz, 0, 1.0, false);
    rk3_step1_corr_forcing_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk1, vrk1, wrk1, 1.0);
    rk3_step1_corr_pressure_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    pflow->u_relax(p,a,pgc,urk1);
    pflow->v_relax(p,a,pgc,vrk1);
    pflow->w_relax(p,a,pgc,wrk1);
    pflow->p_relax(p,a,pgc,a->press);

    rk3_step1_corr_relax_time = pgc->timer() - block_start;

    block_start = pgc->timer();

    #if USE_AMREX
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,urk1,vrk1,wrk1,true);
    pgc->startBatch(p, m_rk1, 0, {{&urk1,gcval_u},{&vrk1,gcval_v},{&wrk1,gcval_w}});
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,urk1,vrk1,wrk1,false);
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

//Step 2
//--------------------------------------------------------

    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);
    irhs(p,a);
    pconvec->start(p,a,urk1,1,urk1,vrk1,wrk1);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk1,urk1,vrk1,wrk1,0.25);

    // ULOOP
    // urk2(i,j,k) = 0.75*a->u(i,j,k) + 0.25*udiff(i,j,k)
    //             + 0.25*p->dt*CPOR1*a->F(i,j,k);
    FIELDLOOP(
        urk2,
        FIELD_CONST_MEMBER(a, u); FIELD_CONST(udiff); FIELD_CONST_MEMBER(a, F); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor1 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i+1,j,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        urk2(i,j,k) = 0.75*member_u(i,j,k) + 0.25*udiff(i,j,k) + 0.25*p->dt*cpor1*member_F(i,j,k);
    )

    p->utime+=pgc->timer()-starttime;
    rk3_step2_u_time = p->utime - rk3_step1_u_time;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    jrhs(p,a);
    pconvec->start(p,a,vrk1,2,urk1,vrk1,wrk1);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk1,urk1,vrk1,wrk1,0.25);

    // VLOOP
    // vrk2(i,j,k) = 0.75*a->v(i,j,k) + 0.25*vdiff(i,j,k)
    //             + 0.25*p->dt*CPOR2*a->G(i,j,k);
    FIELDLOOP(
        vrk2,
        FIELD_CONST_MEMBER(a, v); FIELD_CONST(vdiff); FIELD_CONST_MEMBER(a, G); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor2 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j+1,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        vrk2(i,j,k) = 0.75*member_v(i,j,k) + 0.25*vdiff(i,j,k) + 0.25*p->dt*cpor2*member_G(i,j,k);
    )

    p->vtime+=pgc->timer()-starttime;
    rk3_step2_v_time = p->vtime - rk3_step1_v_time;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    krhs(p,a);
    pconvec->start(p,a,wrk1,3,urk1,vrk1,wrk1);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk1,urk1,vrk1,wrk1,0.25);

    // WLOOP
    // wrk2(i,j,k) = 0.75*a->w(i,j,k) + 0.25*wdiff(i,j,k)
    //             + 0.25*p->dt*CPOR3*a->H(i,j,k);
    FIELDLOOP(
        wrk2,
        FIELD_CONST_MEMBER(a, w); FIELD_CONST(wdiff); FIELD_CONST_MEMBER(a, H); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor3 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j,k+1) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        wrk2(i,j,k) = 0.75*member_w(i,j,k) + 0.25*wdiff(i,j,k) + 0.25*p->dt*cpor3*member_H(i,j,k);
    )

    p->wtime+=pgc->timer()-starttime;
    rk3_step2_w_time = p->wtime - rk3_step1_w_time;

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           urk2, vrk2, wrk2, fx, fy, fz, 1, 0.25, false);

    pflow->pressure_io(p,a,pgc);
    velpat_probe(p,a,wrk2,"s2-pre-proj");
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow, urk2, vrk2, wrk2, 0.25);
    velpat_probe(p,a,wrk2,"s2-post-proj");

    pflow->u_relax(p,a,pgc,urk2);
    pflow->v_relax(p,a,pgc,vrk2);
    pflow->w_relax(p,a,pgc,wrk2);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,urk2,vrk2,wrk2,true);
    pgc->startBatch(p, m_rk2, 0, {{&urk2,gcval_u},{&vrk2,gcval_v},{&wrk2,gcval_w}});
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,urk2,vrk2,wrk2,false);
    #else
    pgc->start1(p,urk2,gcval_u);
    pgc->start2(p,vrk2,gcval_v);
    pgc->start3(p,wrk2,gcval_w);
    #endif
    velpat_probe(p,a,wrk2,"s2-post-gc");
    rk3_step2_corr_time = pgc->timer() - block_start;

//Step 3
//--------------------------------------------------------

    // --- Gated warm-started pressure corrector (env REEF_PCORR_MAX, default 0 = off) ---
    // The stage-3 predictor applies the stage-2 pressure -- one substage stale relative to the
    // current density. At a moving free surface that press/density inconsistency leaves a
    // rotational velocity the single projection cannot remove (AMR_CF_PROJECTION_CONSISTENCY_
    // RECORD.md, Open item #4). Re-run the stage-3 predictor+projection from u^n with the
    // just-solved (refreshed) a->press so grad(press)/roface balances gravity; iterate until the
    // rest velocity falls below REEF_PCORR_TOL (or REEF_PCORR_MAX passes). Warm-started: a->press
    // carries over, so each extra Poisson solve is a small, fast correction. urk1/vrk1/wrk1 are
    // free after step 2 and reused here to hold u^n. Intended for the at-rest / free-surface
    // hydrostatic problem; 6DOF/FSI re-forcing per pass is untested.
    const int    pcorr_max = (std::getenv("REEF_PCORR_MAX") ? std::atoi(std::getenv("REEF_PCORR_MAX")) : 0);
    const double pcorr_tol = (std::getenv("REEF_PCORR_TOL") ? std::atof(std::getenv("REEF_PCORR_TOL")) : 0.0);

    if(pcorr_max>0)
    {
        ULOOP urk1(i,j,k)=a->u(i,j,k);   // snapshot u^n (a->u/v/w still hold it here)
        VLOOP vrk1(i,j,k)=a->v(i,j,k);
        WLOOP wrk1(i,j,k)=a->w(i,j,k);
    }

    for(int pcorr=0; pcorr<=pcorr_max; ++pcorr)
    {
        if(pcorr>0)
        {
            // restart the predictor from u^n; refresh halos so convection sees valid ghosts
            ULOOP a->u(i,j,k)=urk1(i,j,k);
            VLOOP a->v(i,j,k)=vrk1(i,j,k);
            WLOOP a->w(i,j,k)=wrk1(i,j,k);
            #if USE_AMREX
            pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
            #else
            pgc->start1(p,a->u,gcval_u);
            pgc->start2(p,a->v,gcval_v);
            pgc->start3(p,a->w,gcval_w);
            #endif
        }

    // U
    starttime=pgc->timer();

    pturb->isource(p,a);
    pflow->isource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a,a->eta,a->eta_n);
    irhs(p,a);
    pconvec->start(p,a,urk2,1,urk2,vrk2,wrk2);
    pdiff->diff_u(p,a,pgc,psolv,udiff,urk2,urk2,vrk2,wrk2,2.0/3.0);

    // ULOOP
    // a->u(i,j,k) = (1.0/3.0)*a->u(i,j,k) + (2.0/3.0)*udiff(i,j,k)
    //             + (2.0/3.0)*p->dt*CPOR1*a->F(i,j,k);
    FIELDLOOP_MEMBER(a,u,
        FIELD_CONST(udiff); FIELD_CONST_MEMBER(a, F); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor1 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i+1,j,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        u(i,j,k) = (1.0/3.0)*u(i,j,k) + (2.0/3.0)*udiff(i,j,k) + (2.0/3.0)*p->dt*cpor1*member_F(i,j,k);
    )

    p->utime+=pgc->timer()-starttime;
    rk3_step3_u_time = p->utime - rk3_step1_u_time - rk3_step2_u_time;

    // V
    starttime=pgc->timer();

    pturb->jsource(p,a);
    pflow->jsource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a,a->eta,a->eta_n);
    jrhs(p,a);
    pconvec->start(p,a,vrk2,2,urk2,vrk2,wrk2);
    pdiff->diff_v(p,a,pgc,psolv,vdiff,vrk2,urk2,vrk2,wrk2,2.0/3.0);

    // VLOOP
    // a->v(i,j,k) = (1.0/3.0)*a->v(i,j,k) + (2.0/3.0)*vdiff(i,j,k)
    //             + (2.0/3.0)*p->dt*CPOR2*a->G(i,j,k);
    FIELDLOOP_MEMBER(a,v,
        FIELD_CONST(vdiff); FIELD_CONST_MEMBER(a, G); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor2 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j+1,k) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        v(i,j,k) = (1.0/3.0)*v(i,j,k) + (2.0/3.0)*vdiff(i,j,k) + (2.0/3.0)*p->dt*cpor2*member_G(i,j,k);
    )

    p->vtime+=pgc->timer()-starttime;
    rk3_step3_v_time = p->vtime - rk3_step1_v_time - rk3_step2_v_time;

    // W
    starttime=pgc->timer();

    pturb->ksource(p,a);
    pflow->ksource(p,a,pgc,pvrans);
    bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a,a->eta,a->eta_n);
    krhs(p,a);
    pconvec->start(p,a,wrk2,3,urk2,vrk2,wrk2);
    pdiff->diff_w(p,a,pgc,psolv,wdiff,wrk2,urk2,vrk2,wrk2,2.0/3.0);

    // WLOOP
    // a->w(i,j,k) = (1.0/3.0)*a->w(i,j,k) + (2.0/3.0)*wdiff(i,j,k)
    //             + (2.0/3.0)*p->dt*CPOR3*a->H(i,j,k);
    FIELDLOOP_MEMBER(a,w,
        FIELD_CONST(wdiff); FIELD_CONST_MEMBER(a, H); FIELD_CONST_MEMBER(a, porosity);,
        const double cpor3 = (p->B260==0.0 ? 1.0 : 1.0/(1.0+(p->B260*(0.5*(member_porosity(i,j,k+1) + member_porosity(i,j,k))<1.0?1.0:0.0))));
        w(i,j,k) = (1.0/3.0)*w(i,j,k) + (2.0/3.0)*wdiff(i,j,k) + (2.0/3.0)*p->dt*cpor3*member_H(i,j,k);
    )

    p->wtime+=pgc->timer()-starttime;
    rk3_step3_w_time = p->wtime - rk3_step1_w_time - rk3_step2_w_time;

    velpat_probe(p,a,a->w,"s3-pred");   // raw stage-3 predictor output (before forcing/projection)

    block_start = pgc->timer();
    momentum_forcing_start(p, a, pgc, p6dof, pfsi,
                           a->u, a->v, a->w, fx, fy, fz, 2, 2.0/3.0, true);

    pflow->pressure_io(p,a,pgc);
    velpat_probe(p,a,a->w,"s3-pre-proj");
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,2.0/3.0);
    velpat_probe(p,a,a->w,"s3-post-proj");

        if(pcorr_max>0)
        {
            double pcorr_um=0.0;
            ULOOP pcorr_um=std::max(pcorr_um,fabs(a->u(i,j,k)));
            pcorr_um=pgc->globalmax(pcorr_um);
            if(p->mpirank==0)
                std::cout<<"  [pcorr] pass "<<pcorr<<"  umax="<<pcorr_um<<std::endl;
            if(pcorr>0 && pcorr_um<pcorr_tol) break;   // at-rest convergence
        }
    } // --- end gated warm-started pressure corrector loop ---

    pflow->u_relax(p,a,pgc,a->u);
    pflow->v_relax(p,a,pgc,a->v);
    pflow->w_relax(p,a,pgc,a->w);
    pflow->p_relax(p,a,pgc,a->press);

    #if USE_AMREX
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,a->u,a->v,a->w,true);
    pgc->startBatch(p, a->m_mf, 0, {{&a->u,gcval_u},{&a->v,gcval_v},{&a->w,gcval_w}});
    if(cf_lowface_restore) ppoissonsolv->cf_lowface_save_restore(p,a->u,a->v,a->w,false);
    #else
    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
    #endif

    rk3_step3_corr_time = pgc->timer() - block_start;

    const double rk3_total_time = pgc->timer() - rk3_total_start;
    const double rk3_measured_time = rk3_setup_time
                                   + rk3_step1_u_time + rk3_step1_v_time + rk3_step1_w_time + rk3_step1_corr_time
                                   + rk3_step2_u_time + rk3_step2_v_time + rk3_step2_w_time + rk3_step2_corr_time
                                   + rk3_step3_u_time + rk3_step3_v_time + rk3_step3_w_time + rk3_step3_corr_time;
    const double rk3_other_time = rk3_total_time - rk3_measured_time;

    if(p->mpirank==0 && std::getenv("REEF_timing"))
    {
        const double denom = (rk3_total_time > 0.0) ? rk3_total_time : 1.0;
        std::cout<<"momentum RK3 runtime breakdown (s | %total):"<<std::endl;
        std::cout<<std::setprecision(6)
                 <<"  setup: "<<rk3_setup_time<<" | "<<(100.0*rk3_setup_time/denom)<<std::endl
                 <<"  step1 U: "<<rk3_step1_u_time<<" | "<<(100.0*rk3_step1_u_time/denom)<<std::endl
                 <<"  step1 V: "<<rk3_step1_v_time<<" | "<<(100.0*rk3_step1_v_time/denom)<<std::endl
                 <<"  step1 W: "<<rk3_step1_w_time<<" | "<<(100.0*rk3_step1_w_time/denom)<<std::endl
                 <<"  step1 corr: "<<rk3_step1_corr_time<<" | "<<(100.0*rk3_step1_corr_time/denom)<<std::endl
                 <<"    step1 U source: "<<rk3_step1_u_source_time<<" | "<<(100.0*rk3_step1_u_source_time/denom)<<std::endl
                 <<"    step1 U bcmom: "<<rk3_step1_u_bcmom_time<<" | "<<(100.0*rk3_step1_u_bcmom_time/denom)<<std::endl
                 <<"    step1 U pgrad: "<<rk3_step1_u_pgrad_time<<" | "<<(100.0*rk3_step1_u_pgrad_time/denom)<<std::endl
                 <<"    step1 U rhs: "<<rk3_step1_u_rhs_time<<" | "<<(100.0*rk3_step1_u_rhs_time/denom)<<std::endl
                 <<"    step1 U convec: "<<rk3_step1_u_convec_time<<" | "<<(100.0*rk3_step1_u_convec_time/denom)<<std::endl
                 <<"    step1 U diff: "<<rk3_step1_u_diff_time<<" | "<<(100.0*rk3_step1_u_diff_time/denom)<<std::endl
                 <<"    step1 U update: "<<rk3_step1_u_update_time<<" | "<<(100.0*rk3_step1_u_update_time/denom)<<std::endl
                 <<"    step1 V source: "<<rk3_step1_v_source_time<<" | "<<(100.0*rk3_step1_v_source_time/denom)<<std::endl
                 <<"    step1 V bcmom: "<<rk3_step1_v_bcmom_time<<" | "<<(100.0*rk3_step1_v_bcmom_time/denom)<<std::endl
                 <<"    step1 V pgrad: "<<rk3_step1_v_pgrad_time<<" | "<<(100.0*rk3_step1_v_pgrad_time/denom)<<std::endl
                 <<"    step1 V rhs: "<<rk3_step1_v_rhs_time<<" | "<<(100.0*rk3_step1_v_rhs_time/denom)<<std::endl
                 <<"    step1 V convec: "<<rk3_step1_v_convec_time<<" | "<<(100.0*rk3_step1_v_convec_time/denom)<<std::endl
                 <<"    step1 V diff: "<<rk3_step1_v_diff_time<<" | "<<(100.0*rk3_step1_v_diff_time/denom)<<std::endl
                 <<"    step1 V update: "<<rk3_step1_v_update_time<<" | "<<(100.0*rk3_step1_v_update_time/denom)<<std::endl
                 <<"    step1 W source: "<<rk3_step1_w_source_time<<" | "<<(100.0*rk3_step1_w_source_time/denom)<<std::endl
                 <<"    step1 W bcmom: "<<rk3_step1_w_bcmom_time<<" | "<<(100.0*rk3_step1_w_bcmom_time/denom)<<std::endl
                 <<"    step1 W pgrad: "<<rk3_step1_w_pgrad_time<<" | "<<(100.0*rk3_step1_w_pgrad_time/denom)<<std::endl
                 <<"    step1 W rhs: "<<rk3_step1_w_rhs_time<<" | "<<(100.0*rk3_step1_w_rhs_time/denom)<<std::endl
                 <<"    step1 W convec: "<<rk3_step1_w_convec_time<<" | "<<(100.0*rk3_step1_w_convec_time/denom)<<std::endl
                 <<"    step1 W diff: "<<rk3_step1_w_diff_time<<" | "<<(100.0*rk3_step1_w_diff_time/denom)<<std::endl
                 <<"    step1 W update: "<<rk3_step1_w_update_time<<" | "<<(100.0*rk3_step1_w_update_time/denom)<<std::endl
                 <<"    step1 corr forcing: "<<rk3_step1_corr_forcing_time<<" | "<<(100.0*rk3_step1_corr_forcing_time/denom)<<std::endl
                 <<"    step1 corr pressure: "<<rk3_step1_corr_pressure_time<<" | "<<(100.0*rk3_step1_corr_pressure_time/denom)<<std::endl
                 <<"    step1 corr relax: "<<rk3_step1_corr_relax_time<<" | "<<(100.0*rk3_step1_corr_relax_time/denom)<<std::endl
                 <<"    step1 corr gc: "<<rk3_step1_corr_gc_time<<" | "<<(100.0*rk3_step1_corr_gc_time/denom)<<std::endl
                 <<"  step2 U: "<<rk3_step2_u_time<<" | "<<(100.0*rk3_step2_u_time/denom)<<std::endl
                 <<"  step2 V: "<<rk3_step2_v_time<<" | "<<(100.0*rk3_step2_v_time/denom)<<std::endl
                 <<"  step2 W: "<<rk3_step2_w_time<<" | "<<(100.0*rk3_step2_w_time/denom)<<std::endl
                 <<"  step2 corr: "<<rk3_step2_corr_time<<" | "<<(100.0*rk3_step2_corr_time/denom)<<std::endl
                 <<"  step3 U: "<<rk3_step3_u_time<<" | "<<(100.0*rk3_step3_u_time/denom)<<std::endl
                 <<"  step3 V: "<<rk3_step3_v_time<<" | "<<(100.0*rk3_step3_v_time/denom)<<std::endl
                 <<"  step3 W: "<<rk3_step3_w_time<<" | "<<(100.0*rk3_step3_w_time/denom)<<std::endl
                 <<"  step3 corr: "<<rk3_step3_corr_time<<" | "<<(100.0*rk3_step3_corr_time/denom)<<std::endl
                 <<"  other: "<<rk3_other_time<<" | "<<(100.0*rk3_other_time/denom)<<std::endl
                 <<"  total: "<<rk3_total_time<<std::endl;
    }
}

void momentum_RK3::irhs(lexer *p, fdm *a)
{
    n=0;
    const bool relPressure = p->Y9;
    ULOOP
    {
        const double gi = relPressure ? 0.0 : a->gi;
        a->maxF=std::max(fabs(a->rhsvec.V[n] + gi),a->maxF);
        a->F(i,j,k) += (a->rhsvec.V[n] + gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;

        a->rhsvec.V[n]=0.0;
        a->Fext(i,j,k)=0.0;
        ++n;
    }
}

void momentum_RK3::jrhs(lexer *p, fdm *a)
{
    n=0;
    const bool relPressure = p->Y9;
    VLOOP
    {
        const double gj = relPressure ? 0.0 : a->gj;
        a->maxG=std::max(fabs(a->rhsvec.V[n] + gj),a->maxG);
        a->G(i,j,k) += (a->rhsvec.V[n] + gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;

        a->rhsvec.V[n] = 0.0;
        a->Gext(i,j,k) = 0.0;
        ++n;
    }
}

void momentum_RK3::krhs(lexer *p, fdm *a)
{
    n=0;
    const bool relPressure = p->Y9;
    const bool krhs_probe = (std::getenv("REEF_KRHS_PROBE")!=nullptr);
    WLOOP
    {
        const double gk = relPressure ? 0.0 : a->gk;
        a->maxH=std::max(fabs(a->rhsvec.V[n] + gk),a->maxH);
#if USE_AMREX
        if(krhs_probe && p->level==1 && i+p->amr_tile_lo.x==38 && j+p->amr_tile_lo.y==34
           && k+p->amr_tile_lo.z>=20 && k+p->amr_tile_lo.z<=23)
            std::cout<<"  [krhsprobe] rank"<<p->mpirank<<" ("<<i+p->amr_tile_lo.x<<","
                     <<j+p->amr_tile_lo.y<<","<<k+p->amr_tile_lo.z<<") n="<<n
                     <<" H_before="<<a->H(i,j,k)<<" rhsvec="<<a->rhsvec.V[n]
                     <<" gk="<<gk<<" W29z="<<p->W29_z<<" Hext="<<a->Hext(i,j,k)<<std::endl;
#endif
        a->H(i,j,k) += (a->rhsvec.V[n] + gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;

        a->rhsvec.V[n] = 0.0;
        a->Hext(i,j,k) = 0.0;
        ++n;
    }
}
