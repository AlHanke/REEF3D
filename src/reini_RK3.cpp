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

#include"reini_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"

reini_RK3::reini_RK3(lexer* p, int type) : frk1(p), frk2(p), dt(p)
{
    if(p->F50==1)
        gcval_phi=51;
    else if(p->F50==2)
        gcval_phi=52;
    else if(p->F50==3)
        gcval_phi=53;
    else if(p->F50==4)
        gcval_phi=54;

    gcval_iniphi=50;

    if((p->F61>1.0e-20 || p->F60>1.0e-20) && p->F50==1)
        gcval_iniphi=51;
    else if((p->F62>1.0e-20 || p->F60>1.0e-20) && p->F50==2)
        gcval_iniphi=52;
    else if(((p->F61>1.0e-20 && p->F62>1.0e-20) || p->F60>1.0e-20) && p->F50==3)
        gcval_iniphi=53;
    if(type==41)
        gcval_iniphi=50;


    if(p->F46==2)
        ppicard = new picard_f(p);
    else if(p->F46==3)
        ppicard = new picard_lsm(p);
    else
        ppicard = new picard_void(p);

    prdisc = new reinidisc_f(p);

    time_preproc(p);
}

reini_RK3::~reini_RK3()
{
    delete ppicard;
    delete prdisc;
}

void reini_RK3::start(fdm *a, lexer *p, field &f, ghostcell *pgc, ioflow* pflow)
{
    double picardtime, flowtime, rdisctime, calctime, gctime, temptime, picardtime2;
    starttime=pgc->timer();

    if(p->count==0)
        gcval = gcval_iniphi;
    else if(p->count>0)
        gcval = gcval_phi;

    ppicard->volcalc(p,a,pgc,f);

    picardtime=pgc->timer();

    if(p->count==0)
    {
        if(p->mpirank==0)
            cout<<"initializing level set..."<<endl<<endl;
        reiniter=2*int(p->maxlength/(p->F43*p->DXM));
        pgc->start4(p,f,gcval_iniphi);
    }
    else if(p->count>0)
        step(p);

    pflow->fsfrkin(p,a,pgc,frk1);
    pflow->fsfrkin(p,a,pgc,frk2);
    pflow->fsfrkout(p,a,pgc,frk1);
    pflow->fsfrkout(p,a,pgc,frk2);

    flowtime=pgc->timer();

    // Convergence threshold: stop when max|sign(phi0)*(|grad phi|-1)| < tol.
    // The RHS L computed by prdisc is dimensionless; 1e-4 means the level set
    // deviates from a signed-distance function by less than 0.01% per iteration.
    static constexpr double reini_conv_tol = 1.0e-4;

    int iters_done = reiniter;
    for(int q=0;q<reiniter;++q)
    {
        // Step 1
        prdisc->start(p,a,pgc,f,a->L,4);

        if(q==0) rdisctime=pgc->timer();

        FIELDLOOP(
            frk1,
            FIELD_CONST(f); FIELD_CONST(dt); FIELD_CONST_MEMBER(a, L),
            frk1(i,j,k) = f(i,j,k) + dt(i,j,k)*member_L(i,j,k);
        )

        if(q==0) calctime=pgc->timer();

        pgc->start4(p,frk1,gcval);

        if(q==0) gctime=pgc->timer();

        // Step 2
        prdisc->start(p,a,pgc,frk1,a->L,4);

        FIELDLOOP(
            frk2,
            FIELD_CONST(f); FIELD_CONST(frk1); FIELD_CONST(dt); FIELD_CONST_MEMBER(a, L),
            frk2(i,j,k) = 0.75*f(i,j,k) + 0.25*frk1(i,j,k) + 0.25*dt(i,j,k)*member_L(i,j,k);
        )

        pgc->start4(p,frk2,gcval);

        // Step 3
        prdisc->start(p,a,pgc,frk2,a->L,4);

        // Check convergence on the step-3 RHS before applying the update.
        // If L_inf is already below tolerance on iteration q>0 the update would
        // be negligible, so we can skip it and exit early.
        double L_inf = 0.0;
        #if USE_AMREX
        for (int lev = 0; lev < p->nlevs; ++lev)
            L_inf = std::max(L_inf, a->L.GetMultiFab(lev).norm0(0, 0, true));
        #else
        LOOP
            L_inf = MAX(L_inf, fabs(a->L(i,j,k)));
        #endif
        L_inf = pgc->globalmax(L_inf);

        FIELDLOOP(
            f,
            FIELD_CONST(frk2); FIELD_CONST(dt); FIELD_CONST_MEMBER(a, L),
            f(i,j,k) = (1.0/3.0)*f(i,j,k) + (2.0/3.0)*frk2(i,j,k) + (2.0/3.0)*dt(i,j,k)*member_L(i,j,k);
        )

        pgc->start4(p,f,gcval);

        if (q > 0 && L_inf < reini_conv_tol)
        {
            iters_done = q + 1;
            break;
        }
    }

    temptime=pgc->timer();

    ppicard->correct_ls(p,a,pgc,f);

    picardtime2=pgc->timer();

    p->reinitime+=pgc->timer()-starttime;
    if(p->mpirank==0 && p->count>0 && p->count%p->P12==0)
    {
        std::cout<<"\nReini RK3\n"
        <<"  Time for Picard iteration: "<<(picardtime-starttime)*1000.0<<" ms\n"
        <<"  Time for flow calculations: "<<(flowtime-picardtime)*1000.0<<" ms\n"
        <<"  Time for reinitialization discretization: "<<(rdisctime-picardtime)*1000.0<<" ms\n"
        <<"  Time for RK3 calculations: "<<(calctime-rdisctime)*1000.0<<" ms\n"
        <<"  Time for ghost cell updates: "<<(gctime-calctime)*1000.0<<" ms\n"
        <<"  Time for correct_ls: "<<(temptime-picardtime2)*1000.0<<" ms\n"
        <<"  Picard iterations: "<<iters_done<<" / "<<reiniter<<"\n"<<std::endl;
    }
}

void reini_RK3::step(lexer* p)
{
    reiniter=p->F44;
}

void reini_RK3::time_preproc(lexer* p)
{
    #if USE_AMREX
    const int max_level = p->nlevs - 1;
    #else
    const int max_level = 0;
    #endif

    for(int lev=max_level; lev>=0; --lev)
    {
        p->level = lev;
        TILE_LOOP
        IJKLOOP
        PCHECK
        {
            if(p->j_dir==0)
                dt(i,j,k) = p->F43*MIN(p->DXP[IP],p->DZP[KP]);
            else if(p->j_dir==1)
                dt(i,j,k) = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
        }
    }
    p->level = 0;
}
