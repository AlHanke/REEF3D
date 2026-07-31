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
#include<cmath>
#include<cstddef>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<cstdio>

#if USE_AMREX
// Sync the coarse velocity under each fine patch with the fine solution: build face
// MultiFabs from the staggered velocity, average_down_faces (fine -> covered coarse faces,
// incl. the C-F interface face), copy the synced coarse faces back. This makes the coarse
// C-F face velocity equal the sum/average of the fine sub-fluxes (reflux) and the covered
// coarse cells equal the fine average. Adapted from solver_mlmg::average_down_velocity but
// operating on the velocity fields actually being projected. Gated on nlevs>1 by the caller.
static void cf_average_down_velocity(lexer* p, field& u, field& v, field& w)
{
    const int jdir = p->j_dir;
    field* vfld[AMREX_SPACEDIM] = {&u, &v, &w};

    for (int lev = p->nlevs-1; lev >= 1; --lev)
    {
        const int clev = lev-1;
        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> ffine, fcrse;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            ffine[dir].define(amrex::convert(p->amrex_box_array[lev],  e), p->amrex_distribution_mapping[lev],  1, 0);
            fcrse[dir].define(amrex::convert(p->amrex_box_array[clev], e), p->amrex_distribution_mapping[clev], 1, 0);

            // staggered cell-centred velocity -> face MultiFab (face f = vel(f-e))
            auto fill_face = [&] (amrex::MultiFab& fmf, int amrlev)
            {
                auto& v_mf = vfld[dir]->GetMultiFab(amrlev);
                for (amrex::MFIter mfi(fmf); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& fbx = mfi.validbox();
                    auto       f = fmf.array(mfi);
                    const auto vv = v_mf.const_array(mfi);
                    amrex::LoopOnCpu(fbx, [&] (int i, int j, int k) noexcept
                    { f(i,j,k) = vv(i-e[0], j-e[1], k-e[2]); });
                }
            };
            fill_face(ffine[dir], lev);
            fill_face(fcrse[dir], clev);   // pre-fill so uncovered coarse faces survive
        }

        const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> ffp{&ffine[0], &ffine[1], &ffine[2]};
        const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>       fcp{&fcrse[0], &fcrse[1], &fcrse[2]};

        amrex::average_down_faces(ffp, fcp, p->ref_vec, p->amrex_geometry[clev]);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if (dir==1 && jdir!=1) continue;
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            auto& v_mf = vfld[dir]->GetMultiFab(clev);
            for (amrex::MFIter mfi(v_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                auto       vv = v_mf.array(mfi);
                const auto f  = fcrse[dir].const_array(mfi);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
                { vv(i,j,k) = f(i+e[0], j+e[1], k+e[2]); });
            }
        }
    }
}
#endif

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

    // Composite predictor divergence (D side): make the PREDICTOR normal velocity single-valued
    // at the C-F interface -- set each covered coarse face = area-average of the fine faces
    // across it, so the coarse cell's divergence (rhs below) uses the summed fine flux, exactly
    // the flux L's coarse row couples to (conservative M_cf*V_c = M_fc*V_f). This makes the
    // discrete divergence D the adjoint of the gradient G the matrix used, so D(1/rho)G = L.
    #if USE_AMREX
    if(p->nlevs > 1)
    {
        cf_average_down_velocity(p,uvel,vvel,wvel);

        // average_down updated COARSE valid cells under the fine patch; refresh the halos so a
        // neighbour rank's rhs reads the single-valued interface velocity (not the stale
        // pre-average value).
        vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);
    }
    #endif

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

    // C-F interface correction (gradient/G side): correct the fine C-F sub-faces with the
    // centre-distance gradient to the real coarse pcorr, then reflux the coarse velocity under
    // each fine patch to the area-summed fine flux. Together with the predictor average-down
    // this makes D(1/rho)G = L at the interface.
    #if USE_AMREX
    if(p->nlevs > 1)
    {
        psolv->cf_velocity_correction(p,a,pgc,pcorr,uvel,vvel,wvel,alpha);
        cf_average_down_velocity(p,uvel,vvel,wvel);
    }
    #endif

    if(std::getenv("REEF_PROJ_CHECK"))
    projection_consistency_check(p,a,pgc,psolv,alpha);

    p->poissoniter=p->solveriter;

    p->poissontime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"piter: "<<p->solveriter<<"  pres: "<<setprecision(3)<<p->final_res<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
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

void pjm_corr::projection_consistency_check(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, double alpha)
{
    // Pure operator-consistency test of the projection:  (L + D (1/rho) G) pcorr  must be 0.
    // We measure it on an ARBITRARY analytic pcorr (NOT the solved one) so the result is free of
    // the local solver residual A*pcorr-b, which is large at the stiff C-F rows even at tol 1e-9.
    // Apply the full velocity correction G to a zero base -> dU = -(dt/rho) G pcorr, take its
    // divergence R(dU), and compare against L*pcorr. A non-zero residual marks a face where the
    // velocity-correction gradient is not the adjoint of the matrix row -> spurious-velocity source.
    // pcorr is reset by rhs() next step, so overwriting it here is safe.
    constexpr double PI = 3.14159265358979;
    pcorr.setVal(0.0,true);
    LOOP
    pcorr(i,j,k) = cos(PI*p->pos_x())*cos(PI*p->pos_z())
                 * (p->j_dir ? cos(PI*p->pos_y()) : 1.0);
    pgc->start4(p,pcorr,40);

    field4 u0(p), v0(p), w0(p), Lp(p);
    u0.setVal(0.0,true);
    v0.setVal(0.0,true);
    w0.setVal(0.0,true);
    Lp.setVal(0.0,true);

    ucorr(p,a,u0,alpha);
    vcorr(p,a,v0,alpha);
    wcorr(p,a,w0,alpha);

    #if USE_AMREX
    if(p->nlevs > 1)
    {
        psolv->cf_velocity_correction(p,a,pgc,pcorr,u0,v0,w0,alpha);
        cf_average_down_velocity(p,u0,v0,w0);
    }
    #endif

    // Same halo fill the production rhs sees, so D matches the assembled L across boxes/levels.
    vel_setup(p,a,pgc,u0,v0,w0,alpha);

    psolv->matvec_into(p,a,pgc,Lp,pcorr);

    // C-F tag: true if a same-level face neighbour lies outside this level's patch but inside the
    // domain -- i.e. the cell sits on the fine side of a coarse-fine interface (not a domain wall).
    auto cf_tag = [&](int lev, int ci, int cj, int ck) -> bool
    {
        #if USE_AMREX
        const amrex::Box dom = p->amrex_geometry[lev].Domain();
        const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
        for(int d=0; d<6; ++d)
        {
            if(d>=2 && d<4 && !p->j_dir) continue;
            const amrex::IntVect niv(ci+off[d][0], cj+off[d][1], ck+off[d][2]);
            if(dom.contains(niv) && !p->amrex_box_array[lev].intersects(amrex::Box(niv,niv)))
                return true;
        }
        #endif
        return false;
    };

    double r_int = 0.0, r_bnd = 0.0;
    double wi_res = 0.0; int wi_lev = -1, wi[3] = {-1,-1,-1}, wi_flag = 0; bool wi_cov = false;
    double wb_res = 0.0; int wb_lev = -1, wb[3] = {-1,-1,-1}, wb_flag = 0;
    LOOP
    {
        const double div = -(u0(i,j,k) - u0(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
                           -(p->j_dir?(v0(i,j,k) - v0(i,j-1,k))/(alpha*p->dt*p->DYN[JP]):0.0)
                           -(w0(i,j,k) - w0(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);

        const double res = std::fabs(Lp(i,j,k) + div);

        const bool boundary =
               p->flag4(i-1,j,k) < 0 || p->flag4(i+1,j,k) < 0
            || p->flag4(i,j,k-1) < 0 || p->flag4(i,j,k+1) < 0
            || (p->j_dir && (p->flag4(i,j-1,k) < 0 || p->flag4(i,j+1,k) < 0));

        if(boundary)
        {
            r_bnd = std::max(r_bnd,res);
            if(res > wb_res)
            {
                wb_res = res; wb_lev = p->level;
                wb[0]=i; wb[1]=j; wb[2]=k; wb_flag = p->flag4(i,j,k);
            }
        }
        else
        {
            r_int = std::max(r_int,res);
            if(res > wi_res)
            {
                wi_res = res; wi_lev = p->level;
                wi[0]=i; wi[1]=j; wi[2]=k; wi_flag = p->flag4(i,j,k);
                // covered: this cell's refined footprint intersects the next finer level's grids.
                #if USE_AMREX
                wi_cov = false;
                if(p->level < p->nlevs-1)
                {
                    amrex::Box foot(amrex::IntVect(i,j,k), amrex::IntVect(i,j,k));
                    foot.refine(p->ref_vec);
                    wi_cov = p->amrex_box_array[p->level+1].intersects(foot);
                }
                #endif
            }
        }
    }

    const double g_int = pgc->globalmax(r_int);
    const double g_bnd = pgc->globalmax(r_bnd);
    if(p->mpirank==0)
    std::cout<<"\n  [projcheck] max|L*pcorr + R(dU)|  interior="<<g_int
             <<"  boundary="<<g_bnd<<std::endl;
    if(wi_lev>=0 && wi_res==g_int)
    std::cout<<"  [projcheck] worst interior res="<<wi_res<<" at lev="<<wi_lev
             <<" ("<<wi[0]<<","<<wi[1]<<","<<wi[2]<<")  flag4="<<wi_flag
             <<"  covered="<<wi_cov<<"  cf="<<cf_tag(wi_lev,wi[0],wi[1],wi[2])<<std::endl;
    if(wb_lev>=0 && wb_res==g_bnd)
    std::cout<<"  [projcheck] worst boundary res="<<wb_res<<" at lev="<<wb_lev
             <<" ("<<wb[0]<<","<<wb[1]<<","<<wb[2]<<")  flag4="<<wb_flag
             <<"  cf="<<cf_tag(wb_lev,wb[0],wb[1],wb[2])<<std::endl;
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
