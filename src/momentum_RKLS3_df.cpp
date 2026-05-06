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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"momentum_RKLS3_df.h"
#include"vrans.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"ediff2.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"6DOF.h"
#include"FSI.h"

momentum_RKLS3_df::momentum_RKLS3_df
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
):bcmom(p),
    #if USE_AMREX
    m_rk(make_mf(p,3,&m_rk)), m_f(make_mf(p,3,&m_f)),
    urk(p,&m_rk,0), fx(p,&m_f,0),
    vrk(p,&m_rk,1), fy(p,&m_f,1),
    wrk(p,&m_rk,2), fz(p,&m_f,2),
    #else
    urk(p),vrk(p),wrk(p),
    fx(p),fy(p),fz(p),
    #endif
    Cu(p),Cv(p),Cw(p)
{
    pconvec=pconvection;
    pdiff=pdiffusion;
    ppress=ppressure;
    ppois=ppoisson;
    pturb=pturbulence;
    psolv=psolver;
    ppoissonsolv=ppoissonsolver;
    pflow=pioflow;

    alpha << 4.0/15.0, 1.0/15.0, 1.0/6.0;
    gamma << 8.0/15.0, 5.0/12.0, 3.0/4.0;
    zeta << 0.0, -17.0/60.0, -5.0/12.0;
}

void momentum_RKLS3_df::starti(lexer* p, fdm* a, ghostcell* pgc, sixdof* p6dof, vrans* pvrans, fsi* pfsi)
{
    // Set inflow
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);

    bool final = false;

    for (int loop=0; loop<3; loop++)
    {
        if (loop==2) final = true;

        pflow->rkinflow(p,a,pgc,urk,vrk,wrk);

    // -------------------
        // U
        double starttime=pgc->timer();

        // Fill F
        pturb->isource(p,a);
        pflow->isource(p,a,pgc,pvrans);
        bcmom_start(a,p,pgc,pturb,a->u,gcval_u);
        ppress->upgrad(p,a,a->eta,a->eta_n);
        irhs(p,a);
        pdiff->diff_u(p,a,pgc,psolv,urk,a->u,a->u,a->v,a->w,2.0*alpha(loop));

        ULOOP
        urk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR1*a->F(i,j,k);

        // Add convection
        ULOOP
        a->F(i,j,k)=0.0;

        pconvec->start(p,a,a->u,1,a->u,a->v,a->w);

        ULOOP
        urk(i,j,k) += gamma(loop)*p->dt*CPOR1*a->F(i,j,k) + zeta(loop)*p->dt*CPOR1*Cu(i,j,k);

        ULOOP
        Cu(i,j,k)=a->F(i,j,k);

        p->utime+=pgc->timer()-starttime;

    // -------------------
        // V
        starttime=pgc->timer();

        // Add source
        pturb->jsource(p,a);
        pflow->jsource(p,a,pgc,pvrans);
        bcmom_start(a,p,pgc,pturb,a->v,gcval_v);
        ppress->vpgrad(p,a,a->eta,a->eta_n);
        jrhs(p,a);
        pdiff->diff_v(p,a,pgc,psolv,vrk,a->v,a->u,a->v,a->w,2.0*alpha(loop));

        VLOOP
        vrk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR2*a->G(i,j,k);

        // Add convection
        VLOOP
        a->G(i,j,k)=0.0;

        pconvec->start(p,a,a->v,2,a->u,a->v,a->w);

        VLOOP
        vrk(i,j,k) += gamma(loop)*p->dt*CPOR2*a->G(i,j,k) + zeta(loop)*p->dt*CPOR2*Cv(i,j,k);

        VLOOP
        Cv(i,j,k)=a->G(i,j,k);

        p->vtime+=pgc->timer()-starttime;

    // -------------------
        // W
        starttime=pgc->timer();

        // Add source
        pturb->ksource(p,a);
        pflow->ksource(p,a,pgc,pvrans);
        bcmom_start(a,p,pgc,pturb,a->w,gcval_w);
        ppress->wpgrad(p,a,a->eta,a->eta_n);
        krhs(p,a);
        pdiff->diff_w(p,a,pgc,psolv,wrk,a->w,a->u,a->v,a->w,2.0*alpha(loop));

        WLOOP
        wrk(i,j,k) += 2.0*alpha(loop)*p->dt*CPOR3*a->H(i,j,k);

        // Add convection
        WLOOP
        a->H(i,j,k)=0.0;

        pconvec->start(p,a,a->w,3,a->u,a->v,a->w);

        WLOOP
        wrk(i,j,k) += gamma(loop)*p->dt*CPOR3*a->H(i,j,k) + zeta(loop)*p->dt*CPOR3*Cw(i,j,k);

        WLOOP
        Cw(i,j,k)=a->H(i,j,k);

        p->wtime+=pgc->timer()-starttime;

        #if USE_AMREX
        pgc->startBatch(p, m_rk, 0, {{&urk,gcval_u},{&vrk,gcval_v},{&wrk,gcval_w}});
        #else
        pgc->start1(p,urk,gcval_u);
        pgc->start2(p,vrk,gcval_v);
        pgc->start3(p,wrk,gcval_w);
        #endif

    // ----------------------------------------------
        starttime=pgc->timer();
        // Forcing
        fx.setVal(0.0, true);
        fy.setVal(0.0, true);
        fz.setVal(0.0, true);

        pgc->solid_forcing(p,a,2.0*alpha(loop),urk,vrk,wrk,fx,fy,fz);

        p6dof->start_cfd(p,a,pgc,loop,urk,vrk,wrk,fx,fy,fz,final);

        pfsi->forcing(p,a,pgc,2.0*alpha(loop),urk,vrk,wrk,fx,fy,fz,final);

        ULOOP
        {
            a->u(i,j,k) = urk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR1*fx(i,j,k);

            if(p->count<10)
            a->maxF = std::max(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), a->maxF);

            p->fbmax = std::max(fabs(2.0*alpha(loop)*CPOR1*fx(i,j,k)), p->fbmax);
        }

        VLOOP
        {
            a->v(i,j,k) = vrk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR2*fy(i,j,k);

            if(p->count<10)
            a->maxG = std::max(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), a->maxG);

            p->fbmax = std::max(fabs(2.0*alpha(loop)*CPOR2*fy(i,j,k)), p->fbmax);
        }

        WLOOP
        {
            a->w(i,j,k) = wrk(i,j,k) + 2.0*alpha(loop)*p->dt*CPOR3*fz(i,j,k);

            if(p->count<10)
            a->maxH = std::max(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), a->maxH);

            p->fbmax = std::max(fabs(2.0*alpha(loop)*CPOR3*fz(i,j,k)), p->fbmax);
        }

        p->fbtime+=pgc->timer()-starttime;

    // ----------------------------------------------

        // Pressure
        pflow->pressure_io(p,a,pgc);
        ppress->start(a,p,ppois,ppoissonsolv,pgc, pflow, a->u, a->v, a->w, 2.0*alpha(loop));

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
    }
}

void momentum_RKLS3_df::irhs(lexer *p, fdm *a)
{
    n=0;

    ULOOP
    {
        a->maxF = std::max(fabs(a->rhsvec.V[n] + a->gi), a->maxF);
        a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;

        a->rhsvec.V[n] = 0.0;
        a->Fext(i,j,k) = 0.0;
        ++n;
    }

}

void momentum_RKLS3_df::jrhs(lexer *p, fdm *a)
{
    n=0;

    VLOOP
    {
        a->maxG = std::max(fabs(a->rhsvec.V[n] + a->gj), a->maxG);
        a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;

        a->rhsvec.V[n] = 0.0;
        a->Gext(i,j,k) = 0.0;
        ++n;
    }
}

void momentum_RKLS3_df::krhs(lexer *p, fdm *a)
{
    n=0;

    WLOOP
    {
        a->maxH = std::max(fabs(a->rhsvec.V[n] + a->gk), a->maxH);
        a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;

        a->rhsvec.V[n] = 0.0;
        a->Hext(i,j,k) = 0.0;
        ++n;
    }
}
