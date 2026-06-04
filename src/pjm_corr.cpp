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

#include"pjm_corr.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_sf.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"
#include"density_pst.h"
#include<algorithm>
#include<cstddef>

pjm_corr::pjm_corr(lexer* p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc) : pcorr(p), pressure_reference(p)
{
    if(p->F80==0 && p->F300==0 && p->W90==0)
    {
        if(p->W30==0 && p->C10==0 && p->H10==0)
        {
            if(p->Q10==0)
            pd = new density_f(p);
            else if(p->Q10==1)
            pd = new density_f(p);
        }

        if(p->H10==0 && p->W30==1)
        pd = new density_comp(p);

        if(p->H10>0 && p->C10==0)
        pd = new density_heat(p,pheat);

        if(p->C10>0 && p->H10==0)
        pd = new density_conc(p,pconc);
    }

    if(p->F80>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
    pd = new density_vof(p);

    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);

    if(p->F300>0)
    pd = new density_rheo(p);
}

pjm_corr::~pjm_corr()
{
    delete pd;
}

void pjm_corr::start(fdm* a, lexer*p, poisson* ppois, solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    double starttime=pgc->timer();

    vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);

    ppois->start(p,a,pcorr);

    psolv->start(p,a,pgc,pcorr,a->rhsvec,5);

    constexpr int gcval_press = 40;
    pgc->start4(p,pcorr,gcval_press);
    presscorr(p,a,uvel,vvel,wvel,pcorr,alpha);
    reference_start(p,a,pgc);
    pgc->start4(p,a->press,gcval_press);

    ucorr(p,a,uvel,alpha);
    vcorr(p,a,vvel,alpha);
    wcorr(p,a,wvel,alpha);

    // --- DIAGNOSTIC: per-level discrete divergence after the velocity
    // correction, split into far-field vs coarse/fine-zone cells. Remove once
    // the C-F interface reconciliation is settled.
    #if USE_AMREX
    for(int lev=0; lev<p->nlevs; ++lev)
    {
        const amrex::Box      dom = p->amrex_geometry[lev].Domain();
        const amrex::BoxArray ba  = p->amr_cell_mf[lev].boxArray();

        // fine side: neighbour lies inside this level's domain but outside its boxes
        auto is_cf=[&](int ii,int jj,int kk){ amrex::IntVect iv(ii,jj,kk); return dom.contains(iv) && !ba.contains(iv); };

        // coarse side: this cell (or a neighbour) is covered by the next finer level
        const bool has_finer = (lev < p->nlevs-1);
        amrex::BoxArray ba_fine;
        if(has_finer) ba_fine = p->amr_cell_mf[lev+1].boxArray();
        const int rr = p->ref_ratio;
        auto covered=[&](int ii,int jj,int kk){ return has_finer && ba_fine.contains(amrex::IntVect(ii*rr,jj*rr,kk*rr)); };

        double dmax_far=0.0, dmax_cf=0.0;
        int fi=-1,fj=-1,fk=-1; // location of the far-field max on this rank
        p->level=lev;
        for(amrex::MFIter mfi(p->amr_cell_mf[lev]); mfi.isValid(); ++mfi)
        {
            amrex::MFIter *saved=p->set_tile_mfi(&mfi);
            KJILOOP
            PCHECK
            {
                const double dv = (uvel(i,j,k)-uvel(i-1,j,k))/p->DXN[IP]
                                + (p->j_dir?(vvel(i,j,k)-vvel(i,j-1,k))/p->DYN[JP]:0.0)
                                + (wvel(i,j,k)-wvel(i,j,k-1))/p->DZN[KP];

                const int gi=p->amr_tile_lo.x+i, gj=p->amr_tile_lo.y+j, gk=p->amr_tile_lo.z+k;

                const bool nearcf =
                       is_cf(gi-1,gj,gk)||is_cf(gi+1,gj,gk)
                     ||is_cf(gi,gj-1,gk)||is_cf(gi,gj+1,gk)
                     ||is_cf(gi,gj,gk-1)||is_cf(gi,gj,gk+1)
                     ||covered(gi,gj,gk)
                     ||covered(gi-1,gj,gk)||covered(gi+1,gj,gk)
                     ||covered(gi,gj-1,gk)||covered(gi,gj+1,gk)
                     ||covered(gi,gj,gk-1)||covered(gi,gj,gk+1);

                if(nearcf) dmax_cf =MAX(dmax_cf, fabs(dv));
                else if(fabs(dv)>dmax_far) { dmax_far=fabs(dv); fi=gi; fj=gj; fk=gk; }
            }
            p->set_tile_mfi(saved?saved:p->default_cell_mfi.get());
        }
        p->level=0;
        p->set_tile_mfi(p->default_cell_mfi.get());

        const double gmax_far=pgc->globalmax(dmax_far);
        const double gmax_cf =pgc->globalmax(dmax_cf);

        if(p->mpirank==0)
        cout<<"\n  [div] level "<<lev<<"  far-field max|div|: "<<gmax_far
            <<"   C-F-zone max|div|: "<<gmax_cf<<endl;

        if(dmax_far==gmax_far && gmax_far>0.0)
        cout<<"        far-field max at global (i,j,k)=("<<fi<<","<<fj<<","<<fk
            <<")  domain hi=("<<dom.bigEnd(0)<<","<<dom.bigEnd(1)<<","<<dom.bigEnd(2)
            <<")  rank "<<p->mpirank<<endl;
    }
    #endif

    p->poissoniter=p->solveriter;

    p->poissontime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void pjm_corr::ucorr(lexer* p, fdm* a, field& uvel, double alpha)
{
    ULOOP
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((pcorr(i+1,j,k)-pcorr(i,j,k))/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
}

void pjm_corr::vcorr(lexer* p, fdm* a, field& vvel, double alpha)
{
    if(p->j_dir==1)
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((pcorr(i,j+1,k)-pcorr(i,j,k))/(p->DYP[JP]*pd->roface(p,a,0,1,0)));
}

void pjm_corr::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{
    WLOOP
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((pcorr(i,j,k+1)-pcorr(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
}

void pjm_corr::presscorr(lexer* p, fdm* a, field& uvel, field& vvel, field& wvel, field& pcorr, double alpha)
{
    LOOP
    a->press(i,j,k) += pcorr(i,j,k);
}

void pjm_corr::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w, double alpha)
{
    double uvel,vvel,wvel;

    std::fill(a->rhsvec.V.begin(),a->rhsvec.V.end(),0.0);

    pcorr.setVal(0.0,true);

    size_t count=0;
    LOOP
    {
        a->rhsvec.V[count] = -(u(i,j,k) - u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
                             -(p->j_dir?(v(i,j,k) - v(i,j-1,k))/(alpha*p->dt*p->DYN[JP]):0.0)
                             -(w(i,j,k) - w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);

        ++count;
    }
}

void pjm_corr::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w, double alpha)
{
    constexpr int gcval_u=7;
    constexpr int gcval_v=8;
    constexpr int gcval_w=9;
    pgc->start1(p,u,gcval_u);
    pgc->start2(p,v,gcval_v);
    pgc->start3(p,w,gcval_w);
}

void pjm_corr::upgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    ULOOP
    a->F(i,j,k) -= PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*pd->roface(p,a,1,0,0));
}

void pjm_corr::vpgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    if(p->j_dir)
    VLOOP
    a->G(i,j,k) -= PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*pd->roface(p,a,0,1,0));
}

void pjm_corr::wpgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    WLOOP
    a->H(i,j,k) -= PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1));
}

void pjm_corr::ini(lexer*p, fdm* a, ghostcell *pgc)
{
    reference_ini(p,a,pgc);
}
