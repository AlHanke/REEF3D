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

// REEF_COVERED_PRESS_RECON: overwrite covered coarse press with the fine average (fine is
// authoritative). The per-column hydrostatic build (ini_hydrostatic) leaves each coarse column a
// slightly different constant offset -- invisible to wpgrad (vertical) but a spurious HORIZONTAL
// gradient upgrad/vpgrad turn into the geyser. Averaging the fine solution onto covered coarse
// cells enforces horizontal consistency AND feeds a clean fine C-F press ghost. Re-enables what
// start4(...,false) disables; average_down writes ONLY covered cells (uncovered have no fine
// coverage), so the uncovered coarse column is untouched. Call BEFORE start4 so the C-F ghost fill
// picks up the corrected covered values.
#if USE_AMREX
static void covered_press_avgdown(lexer* p, field& press)
{
    if(p->nlevs<=1) return;
    for(int lev=p->nlevs-2; lev>=0; --lev)
        amrex::average_down(press.GetMultiFab(lev+1), press.GetMultiFab(lev), 0, 1, p->ref_vec);
}
#endif

// True if the cell (tile-local i,j,k on p->level) is COVERED by the next finer level -- its
// refined footprint intersects level+1's grids. Covered coarse cells are fine-authoritative,
// so their predictor pressure-gradient force is spurious (option REEF_SKIP_COVERED_PGRAD).
static bool cell_is_covered(lexer* p, int ci, int cj, int ck)
{
    if(p->level >= p->nlevs-1) return false;
    const int gi=ci+p->amr_tile_lo.x, gj=cj+p->amr_tile_lo.y, gk=ck+p->amr_tile_lo.z;
    amrex::Box foot(amrex::IntVect(gi,gj,gk), amrex::IntVect(gi,gj,gk));
    foot.refine(p->ref_vec);
    return p->amrex_box_array[p->level+1].intersects(foot);
}
#endif

// Bug #2 probe helper (env REEF_PREDDIV): per-level max|D.u*| of the (predictor) velocity, tagged
// by RK stage (alpha 1.0/0.25/0.667) and a call-site tag. At the worst cell it dumps the pressure
// 7-pt stencil + staggered face velocities + flags, so the injecting term can be identified and
// compared flip vs noflip. No-op unless REEF_PREDDIV is set and nlevs>1. Strip once bug #2 closes.
static void preddiv_probe(lexer* p, fdm* a, ghostcell* pgc,
                          field& uvel, field& vvel, field& wvel, double alpha, const char* tag)
{
#if USE_AMREX
    if(!(std::getenv("REEF_PREDDIV") && p->nlevs>1)) return;
    for(int lev=0; lev<p->nlevs; ++lev)
    {
        const double dx=p->amrex_geometry[lev].CellSize(0);
        const double dy=p->amrex_geometry[lev].CellSize(1);
        const double dz=p->amrex_geometry[lev].CellSize(2);
        auto& u_mf=uvel.GetMultiFab(lev);
        auto& v_mf=vvel.GetMultiFab(lev);
        auto& w_mf=wvel.GetMultiFab(lev);
        auto& pr_mf=a->press.GetMultiFab(lev);
        const auto& f_mf=p->flag4.GetMultiFab(lev);
        double dmax=0.0; int wi=-1,wj=-1,wk=-1;
        double pst[7]={0,0,0,0,0,0,0}, ust[2]={0,0}, wst[2]={0,0}; int fst[3]={0,0,0};
        for(amrex::MFIter mfi(u_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx=mfi.validbox();
            const auto ua=u_mf.const_array(mfi);
            const auto va=v_mf.const_array(mfi);
            const auto wa=w_mf.const_array(mfi);
            const auto pa=pr_mf.const_array(mfi);
            const auto fa=f_mf.const_array(mfi);
            for(int kk=bx.smallEnd(2);kk<=bx.bigEnd(2);++kk)
            for(int jj=bx.smallEnd(1);jj<=bx.bigEnd(1);++jj)
            for(int ii=bx.smallEnd(0);ii<=bx.bigEnd(0);++ii)
            {
                if(fa(ii,jj,kk)<AIR_FLAG) continue;
                double d=(ua(ii,jj,kk)-ua(ii-1,jj,kk))/dx+(wa(ii,jj,kk)-wa(ii,jj,kk-1))/dz;
                if(p->j_dir) d+=(va(ii,jj,kk)-va(ii,jj-1,kk))/dy;
                if(std::fabs(d)>dmax)
                {
                    dmax=std::fabs(d); wi=ii; wj=jj; wk=kk;
                    pst[0]=pa(ii,jj,kk);
                    pst[1]=pa(ii-1,jj,kk); pst[2]=pa(ii+1,jj,kk);
                    pst[3]=pa(ii,jj-1,kk); pst[4]=pa(ii,jj+1,kk);
                    pst[5]=pa(ii,jj,kk-1); pst[6]=pa(ii,jj,kk+1);
                    ust[0]=ua(ii-1,jj,kk); ust[1]=ua(ii,jj,kk);
                    wst[0]=wa(ii,jj,kk-1); wst[1]=wa(ii,jj,kk);
                    fst[0]=fa(ii,jj,kk); fst[1]=fa(ii,jj,kk-1); fst[2]=fa(ii,jj,kk+1);
                }
            }
        }
        const double gd=pgc->globalmax(dmax);
        if(p->mpirank==0)
            std::cout<<"  [preddiv] "<<tag<<" stage(alpha="<<alpha<<") lev="<<lev<<" max|D.u*|="<<gd<<std::endl;
        if(dmax==gd && wi>=0)   // the rank holding the global worst dumps its stencil
            std::cout<<"  [preddiv] "<<tag<<" lev="<<lev<<" worst@("<<wi<<","<<wj<<","<<wk<<") rank="<<p->mpirank
                     <<"  press c="<<pst[0]<<" x-/x+="<<pst[1]<<"/"<<pst[2]
                     <<" y-/y+="<<pst[3]<<"/"<<pst[4]<<" z-/z+="<<pst[5]<<"/"<<pst[6]
                     <<"  w z-/z+="<<wst[0]<<"/"<<wst[1]<<"  u x-/x+="<<ust[0]<<"/"<<ust[1]
                     <<"  flag c/z-/z+="<<fst[0]<<"/"<<fst[1]<<"/"<<fst[2]<<std::endl;
    }
#endif
}

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

void pjm_corr::start(fdm* a, lexer* p, poisson* ppois, solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    double starttime=pgc->timer();

    vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);

    // Bug #2 probe (env REEF_PREDDIV): the predictor divergence D.u* per level, measured BEFORE
    // and AFTER cf_predictor_sync, so we can see whether cf_predictor_sync's fill-from-coarse +
    // reflux is what injects the stage-3 checkerboard (the fine cell shows w z-/z+ = -/+0.872 with a
    // SMOOTH pressure, so it is a velocity mode, not -grad(press)). If "pre-sync" is clean and
    // "post-sync" is huge, cf_predictor_sync is the seeder.
    preddiv_probe(p,a,pgc,uvel,vvel,wvel,alpha,"pre-sync");
    cf_predictor_sync(p,a,pgc,uvel,vvel,wvel,psolv,alpha);
    preddiv_probe(p,a,pgc,uvel,vvel,wvel,alpha,"post-sync");

    rhs(p,a,pgc,uvel,vvel,wvel,alpha);

    ppois->start(p,a,pcorr);

    psolv->start(p,a,pgc,pcorr,a->rhsvec,5);

    // Clear the pressure divergence RHS out of the shared a->rhsvec buffer now that the solve
    // has consumed it. a->rhsvec is filled here by LOOP order (rhs(): -div(u)/(alpha*dt)) but is
    // ALSO read by the momentum predictor (momentum_RK3::irhs/jrhs/krhs) by ULOOP/VLOOP/WLOOP
    // order. Those enumerations only coincide when the cell layout is stable; after an AMR
    // regrid/rank-change they desync, so the next stage's krhs reads a DIFFERENT cell's leftover
    // divergence (~ -div/(alpha*dt), O(100) at a surface/C-F cell) and injects it into a->H --
    // a spurious vertical force that seeds a w checkerboard and blows the run up (N61). The term
    // is ~0 at rest with a matched layout (why single-rank/no-flip is stable), so zeroing it here
    // removes the cross-wired read with no effect on the well-behaved path. (Diffusion, when on,
    // does not rely on this: rhs() above already overwrote any diffusion RHS left in a->rhsvec.)
    std::fill(a->rhsvec.V.begin(), a->rhsvec.V.end(), 0.0);

    // REEF_CF_PROJECTION_GROUP member (5) — selects the matrix-consistent C-F ghost fill.
    // gcv 41: two-point d_cf prolongation so the velocity/pressure-gradient corrector matches
    // the SSAMG C-F stencil. Domain BCs are identical to the pressure Neumann gcv 40 (see
    // field_amrex::FillDomainBoundaryImpl). Canonical note in hypre_ssamg::amr_cf_coefficients.
    #if USE_AMREX
    const int gcval_press = (p->nlevs > 1) ? 41 : 40;
    // Predictor press ghost: gcv 42 = transverse-linear C-F ghost (REEF_CF_GHOST_TRANSVERSE),
    // so the next predictor's -grad(press) has no spurious lateral C-F gradient. pcorr stays 41
    // (matrix-consistent d_cf) so the projection/projcheck is untouched.
    const int gcval_press_pred = (p->nlevs > 1 && p->Y11==1) ? 42 : gcval_press;
    #else
    const int gcval_press = 40;
    const int gcval_press_pred = gcval_press;
    #endif
    pgc->start4(p,pcorr,gcval_press,false);
    #if USE_AMREX
    LEVEL_LOOP
    {
        auto const& test_mf = a->test.GetMultiFab();
        for (amrex::MFIter mfi(a->test.GetMultiFab()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            auto const& cc_arr = a->test.GetMultiFab().array(mfi);
            auto const& p_fc = pcorr.GetMultiFab().const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cc_arr(i,j,k,0) += p_fc(i,j,k);
            });
        }
    }
    #endif
    // amrex::MultiFab::Copy(a->test.GetMultiFab(), pcorr.GetMultiFab(), 0, 0, 1, 0);
    presscorr(p,a);
    reference_start(p,a,pgc);
    #if USE_AMREX
    if(std::getenv("REEF_COVERED_PRESS_RECON")) covered_press_avgdown(p,a->press);
    #endif
    pgc->start4(p,a->press,gcval_press_pred,false);   // no average_down: keep coarse press self-consistent
                                                      // (avg-down breaks the hydrostatic grad at surface)

    velcorr(p,a,pgc,uvel,vvel,wvel,psolv,alpha);

#if USE_AMREX
    // Bug #2 probe (env REEF_PCORR_LEVEL): per-level max|pcorr| after the solve and per-level
    // max|D.u| of the corrected velocity, tagged by RK stage (alpha 1.0/0.25/0.667 = stage 1/2/3).
    // Localises the flip-amplified runaway: which LEVEL's pcorr grows, and at which stage. If the
    // coarse (lev0) pcorr blows up while the fine (lev1) stays bounded, the coarse solve is being
    // polluted by the fine->coarse reflux; if lev1 leads, the fine low-face correction is the seed.
    if(std::getenv("REEF_PCORR_LEVEL") && p->nlevs>1)
    {
        for(int lev=0; lev<p->nlevs; ++lev)
        {
            double pmax=0.0, dmax=0.0; int di=-1,dj=-1,dk=-1;
            const double dx=p->amrex_geometry[lev].CellSize(0);
            const double dy=p->amrex_geometry[lev].CellSize(1);
            const double dz=p->amrex_geometry[lev].CellSize(2);
            auto& pc_mf=pcorr.GetMultiFab(lev);
            auto& u_mf=uvel.GetMultiFab(lev);
            auto& v_mf=vvel.GetMultiFab(lev);
            auto& w_mf=wvel.GetMultiFab(lev);
            const auto& f_mf=p->flag4.GetMultiFab(lev);
            for(amrex::MFIter mfi(pc_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx=mfi.validbox();
                const auto pa=pc_mf.const_array(mfi);
                const auto ua=u_mf.const_array(mfi);
                const auto va=v_mf.const_array(mfi);
                const auto wa=w_mf.const_array(mfi);
                const auto fa=f_mf.const_array(mfi);
                for(int kk=bx.smallEnd(2);kk<=bx.bigEnd(2);++kk)
                for(int jj=bx.smallEnd(1);jj<=bx.bigEnd(1);++jj)
                for(int ii=bx.smallEnd(0);ii<=bx.bigEnd(0);++ii)
                {
                    if(fa(ii,jj,kk)<AIR_FLAG) continue;
                    pmax=std::max(pmax,std::fabs(pa(ii,jj,kk)));
                    double d=(ua(ii,jj,kk)-ua(ii-1,jj,kk))/dx+(wa(ii,jj,kk)-wa(ii,jj,kk-1))/dz;
                    if(p->j_dir) d+=(va(ii,jj,kk)-va(ii,jj-1,kk))/dy;
                    if(std::fabs(d)>dmax){dmax=std::fabs(d);di=ii;dj=jj;dk=kk;}
                }
            }
            const double gp=pgc->globalmax(pmax), gd=pgc->globalmax(dmax);
            if(p->mpirank==0)
                std::cout<<"  [pcorrlev] stage(alpha="<<alpha<<") lev="<<lev
                         <<"  max|pcorr|="<<gp<<"  max|D.u|="<<gd
                         <<"  @("<<di<<","<<dj<<","<<dk<<")"<<std::endl;
        }
    }
#endif

    if(std::getenv("REEF_PROJ_CHECK") || std::getenv("REEF_PROJ_ADJOINT"))
    projection_consistency_check(p,a,pgc,psolv,alpha);

    p->poissoniter=p->solveriter;

    p->poissontime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"piter: "<<p->solveriter<<"  pres: "<<setprecision(3)<<p->final_res<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

// REEF_CF_PROJECTION_GROUP member (6a) — u/v/wcorr apply the fine-spacing gradient
// (pcorr(i+1)-pcorr(i))/(DXP*roface). At the C-F interface the ghost pcorr(i+1) is set by
// FillCoarseFineCellGhost so this DXP-spacing difference equals the SSAMG d_cf flux; the
// high C-F faces themselves are sole-written by cf_velocity_correction. Keep DXP/roface here
// consistent with the group (canonical note in hypre_ssamg::amr_cf_coefficients).
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

void pjm_corr::wcorr(lexer* p, fdm* a, field& wvel, double alpha)
{
    WLOOP
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((pcorr(i,j,k+1)-pcorr(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
}

void pjm_corr::velcorr(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel, solver* psolv, double alpha)
{
    #if USE_AMREX
    if(p->nlevs <= 1)
    #endif
    {
        ucorr(p,a,uvel,alpha);
        vcorr(p,a,vvel,alpha);
        wcorr(p,a,wvel,alpha);
        return;
    }

    #if USE_AMREX
    // Multi-level projection velocity correction is owned end-to-end by the solver
    // (hypre_ssamg::velcorr_amr): the interior gradient is sourced from the assembled operator
    // a->M (exact matrix adjoint, auto-zero at C-F/wall faces -> no cf_masks), the free-surface
    // Dirichlet faces from pd->roface, then the C-F sub-face correction and the covered-coarse
    // reflux. This replaces the former interior ULOOP/VLOOP/WLOOP + cf_masks + cf_velocity_correction
    // + cf_average_down_velocity chain that lived here; keeping it in the solver guarantees the
    // corrector reads the operator (a->M / cf_links / a->Mrow) while it is the pressure operator.
    psolv->velcorr_amr(p,a,pgc,pd,pcorr,uvel,vvel,wvel,alpha);
    #endif
}

void pjm_corr::presscorr(lexer* p, fdm* a)
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
    #if USE_AMREX
    const int gcval_press = (p->nlevs > 1) ? 41 : 40;
    #else
    const int gcval_press = 40;
    #endif
    pgc->start4(p,pcorr,gcval_press);   // match production: gcv 41 does the stagger-correct C-F ghost

    field4 u0(p), v0(p), w0(p), Lp(p);
    u0.setVal(0.0,true);
    v0.setVal(0.0,true);
    w0.setVal(0.0,true);
    Lp.setVal(0.0,true);

    // Mirror production (start) exactly: velcorr already runs the masked interior correction
    // AND the C-F correction + reflux (cf_velocity_correction, cf_average_down_velocity) internally,
    // so it must be called ONCE. (Previously this was followed by a second cf_velocity_correction +
    // cf_average_down_velocity, which double-corrected the C-F faces and reported a spurious ~6x
    // inflated V*div -- a diagnostic artifact, not an operator inconsistency.)
    velcorr(p,a,pgc,u0,v0,w0,psolv,alpha);

    // DIAGNOSTIC (REEF_CF_LOWFACE_RESTORE): the open question is whether preserving the
    // velcorr-corrected fine LOW C-F faces through the start1/2/3 below would make the projection
    // consistent. Save them now, let vel_setup clobber them, then restore -- if projcheck drops to
    // ~1e-16 at the C-F cell the hypothesis holds and this save/restore is the production fix.
    const bool lowface_restore = (std::getenv("REEF_CF_LOWFACE_RESTORE") != nullptr);
    if(lowface_restore) psolv->cf_lowface_save_restore(p,u0,v0,w0,true);

    // Same halo fill the production rhs sees, so D matches the assembled L across boxes/levels.
    vel_setup(p,a,pgc,u0,v0,w0,alpha);

    if(lowface_restore) psolv->cf_lowface_save_restore(p,u0,v0,w0,false);

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

    // Covered coarse cell: its refined footprint intersects the next finer level's grids.
    // These are emitted as identity rows (not part of the composite solution; the fine level
    // is authoritative), so they must be excluded from the projection residual buckets.
    // NOTE: ci,cj,ck arrive as TILE-LOCAL LOOP indices; amrex_box_array is in GLOBAL indices,
    // so add amr_tile_lo (the same offset velcorr applies for cf_masks) before the test --
    // otherwise covered cells in tiles with a nonzero origin are missed and show up as a
    // spurious interior residual (an identity row evaluated as if it were solved).
    auto is_covered = [&](int lev, int ci, int cj, int ck) -> bool
    {
        #if USE_AMREX
        if(lev < p->nlevs-1)
        {
            const int gi=ci+p->amr_tile_lo.x, gj=cj+p->amr_tile_lo.y, gk=ck+p->amr_tile_lo.z;
            amrex::Box foot(amrex::IntVect(gi,gj,gk), amrex::IntVect(gi,gj,gk));
            foot.refine(p->ref_vec);
            return p->amrex_box_array[lev+1].intersects(foot);
        }
        #endif
        return false;
    };

    double r_int = 0.0, r_bnd = 0.0;
    double wi_res = 0.0; int wi_lev = -1, wi[3] = {-1,-1,-1}, wi_flag = 0; bool wi_cov = false;
    double wb_res = 0.0; int wb_lev = -1, wb[3] = {-1,-1,-1}, wb_flag = 0;

    // REEF_PROJ_ADJOINT: prove velcorr_amr is the exact matrix adjoint where it matters most.
    // The interior/boundary split above hides the two cases the C-F reflux can break: a C-F
    // fine cell (all fluid neighbours -> lands in "interior") and a free-surface cell (air
    // neighbour -> lumped with domain walls in "boundary"). Track them as their own buckets so
    // the pass criterion (solver-tol at C-F corners AND the surface band) is read off directly.
    const bool adjoint = (std::getenv("REEF_PROJ_ADJOINT") != nullptr);
    double r_cf = 0.0;   int wc_lev = -1, wc[3] = {-1,-1,-1}, wc_flag = 0;
    double r_surf = 0.0; int ws_lev = -1, ws[3] = {-1,-1,-1}, ws_flag = 0;
    // Per-axis breakdown of the worst C-F cell: which face carries the |Lp+V*div| gap.
    double wc_Lp = 0.0, wc_div = 0.0, wc_pc = 0.0, wc_dax[3] = {0,0,0};
    int    wc_fln[6] = {0,0,0,0,0,0}; bool wc_covn[6] = {false,false,false,false,false,false};

    // Neighbourhood dump at the worst interior cell (REEF_PROJ_PROBE): the per-axis residual
    // split (Lp vs V*div), this cell's pcorr/Lp, and for each of the 6 neighbours its flag4 and
    // whether it is COVERED (refined footprint under the next finer patch). Pinpoints whether the
    // AMR-introduced inconsistency comes from a patch-adjacent face.
    const bool projprobe = (std::getenv("REEF_PROJ_PROBE") != nullptr);
    double wi_div = 0.0, wi_Lp = 0.0, wi_pc = 0.0;
    double wi_dax[3] = {0,0,0};   // per-axis V*div contribution (x,y,z) -- which face carries the residual
    double wi_pcn[6] = {0,0,0,0,0,0};
    int    wi_fln[6] = {0,0,0,0,0,0};
    bool   wi_covn[6] = {false,false,false,false,false,false};
    // Per-face density the CORRECTOR/interior-matrix used (pd->roface) + the phi driving it.
    // Interior faces cancel identically (poisson_pcorr M_f and velcorr both divide by the SAME
    // roface), so a nonzero residual localises to the C-F face; these expose the density and the
    // phi contrast there so it can be compared with the matrix C-F coupling (|L.coeff|, whose
    // implied density is printed by REEF_CF_FACE_PROBE in cf_velocity_correction).
    double wi_rof[6] = {0,0,0,0,0,0};   // roface at (x-,x+,y-,y+,z-,z+)
    double wi_phis = 0.0, wi_phin[6] = {0,0,0,0,0,0};
    auto nb_covered = [&](int lev, int ci, int cj, int ck) -> bool
    {
        #if USE_AMREX
        if(lev < p->nlevs-1)
        {
            const int gi=ci+p->amr_tile_lo.x, gj=cj+p->amr_tile_lo.y, gk=ck+p->amr_tile_lo.z;
            amrex::Box foot(amrex::IntVect(gi,gj,gk), amrex::IntVect(gi,gj,gk));
            foot.refine(p->ref_vec);
            return p->amrex_box_array[lev+1].intersects(foot);
        }
        #endif
        return false;
    };

    LOOP
    {
        if(is_covered(p->level,i,j,k)) continue;   // identity row, not part of the solution

        const double div = -(u0(i,j,k) - u0(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
                           -(p->j_dir?(v0(i,j,k) - v0(i,j-1,k))/(alpha*p->dt*p->DYN[JP]):0.0)
                           -(w0(i,j,k) - w0(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);

        // matvec_into applies the assembled, volume-weighted operator, so Lp = V_lev*L_phys*pcorr.
        // The divergence term is physical (unweighted), so weight it by the same V_lev to compare
        // like with like: a consistent projection gives V_lev*(L_phys*pcorr + div) = 0.
        double V_lev = 1.0;
        #if USE_AMREX
        V_lev = p->amrex_geometry[p->level].CellSize(0)
              * (p->j_dir ? p->amrex_geometry[p->level].CellSize(1) : 1.0)
              * p->amrex_geometry[p->level].CellSize(2);
        #endif

        const double res = std::fabs(Lp(i,j,k) + V_lev*div);

        const bool boundary =
               p->flag4(i-1,j,k) < 0 || p->flag4(i+1,j,k) < 0
            || p->flag4(i,j,k-1) < 0 || p->flag4(i,j,k+1) < 0
            || (p->j_dir && (p->flag4(i,j-1,k) < 0 || p->flag4(i,j+1,k) < 0));

        if(adjoint)
        {
            // Global indices for cf_tag / amrex_box_array (both live in global space).
            #if USE_AMREX
            const int gi=i+p->amr_tile_lo.x, gj=j+p->amr_tile_lo.y, gk=k+p->amr_tile_lo.z;
            #else
            const int gi=i+p->origin_i, gj=j+p->origin_j, gk=k+p->origin_k;
            #endif

            // C-F bucket: either the FINE side (cf_tag: a same-level neighbour lies outside the
            // patch) OR the COARSE side (an uncovered cell with a covered/refined face-neighbour --
            // cf_tag cannot see this, so it otherwise hides in the interior bucket). Both are the
            // reflux-critical faces velcorr_amr (M==0 at the seam) + cf_velocity_correction + the
            // reflux must reproduce, or the interface injects divergence.
            bool coarse_cf = false;
            #if USE_AMREX
            if(p->level < p->nlevs-1)
                coarse_cf = nb_covered(p->level,i-1,j,k) || nb_covered(p->level,i+1,j,k)
                         || nb_covered(p->level,i,j,k-1) || nb_covered(p->level,i,j,k+1)
                         || (p->j_dir && (nb_covered(p->level,i,j-1,k) || nb_covered(p->level,i,j+1,k)));
            #endif
            if((cf_tag(p->level,gi,gj,gk) || coarse_cf) && res > r_cf)
            {
                r_cf = res; wc_lev = p->level; wc[0]=gi; wc[1]=gj; wc[2]=gk; wc_flag = p->flag4(i,j,k);
                // Per-axis V*div split: the residual axis is the C-F (or refluxed) face.
                wc_Lp = Lp(i,j,k); wc_div = V_lev*div; wc_pc = pcorr(i,j,k);
                wc_dax[0] = -V_lev*(u0(i,j,k) - u0(i-1,j,k))/(alpha*p->dt*p->DXN[IP]);
                wc_dax[1] = p->j_dir ? -V_lev*(v0(i,j,k) - v0(i,j-1,k))/(alpha*p->dt*p->DYN[JP]) : 0.0;
                wc_dax[2] = -V_lev*(w0(i,j,k) - w0(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                wc_fln[0]=p->flag4(i-1,j,k); wc_fln[1]=p->flag4(i+1,j,k);
                wc_fln[2]=p->flag4(i,j-1,k); wc_fln[3]=p->flag4(i,j+1,k);
                wc_fln[4]=p->flag4(i,j,k-1); wc_fln[5]=p->flag4(i,j,k+1);
                #if USE_AMREX
                wc_covn[0]=nb_covered(p->level,i-1,j,k); wc_covn[1]=nb_covered(p->level,i+1,j,k);
                wc_covn[2]=nb_covered(p->level,i,j-1,k); wc_covn[3]=nb_covered(p->level,i,j+1,k);
                wc_covn[4]=nb_covered(p->level,i,j,k-1); wc_covn[5]=nb_covered(p->level,i,j,k+1);
                #endif
            }

            // Surface bucket: an AIR_FLAG neighbour -> the Dirichlet face velcorr_amr corrects
            // with the pd->roface fall-back. Separated from domain walls so a surface-only
            // residual is not masked by the (start1/2/3-owned) wall faces in r_bnd.
            const bool surf =
                   p->flag4(i-1,j,k)==AIR_FLAG || p->flag4(i+1,j,k)==AIR_FLAG
                || p->flag4(i,j,k-1)==AIR_FLAG || p->flag4(i,j,k+1)==AIR_FLAG
                || (p->j_dir && (p->flag4(i,j-1,k)==AIR_FLAG || p->flag4(i,j+1,k)==AIR_FLAG));
            if(surf && res > r_surf)
            {
                r_surf = res; ws_lev = p->level; ws[0]=gi; ws[1]=gj; ws[2]=gk; ws_flag = p->flag4(i,j,k);
            }
        }

        if(boundary)
        {
            r_bnd = std::max(r_bnd,res);
            if(res > wb_res)
            {
                wb_res = res; wb_lev = p->level;
                // store GLOBAL indices (cf_tag / the print test amrex_box_array in global space)
                #if USE_AMREX
                wb[0]=i+p->amr_tile_lo.x; wb[1]=j+p->amr_tile_lo.y; wb[2]=k+p->amr_tile_lo.z;
                #else
                wb[0]=i+p->origin_i; wb[1]=j+p->origin_j; wb[2]=k+p->origin_k;
                #endif
                wb_flag = p->flag4(i,j,k);
            }
        }
        else
        {
            r_int = std::max(r_int,res);
            if(res > wi_res)
            {
                wi_res = res; wi_lev = p->level;
                // store GLOBAL indices (cf_tag / the print test amrex_box_array in global space)
                #if USE_AMREX
                wi[0]=i+p->amr_tile_lo.x; wi[1]=j+p->amr_tile_lo.y; wi[2]=k+p->amr_tile_lo.z;
                #else
                wi[0]=i+p->origin_i; wi[1]=j+p->origin_j; wi[2]=k+p->origin_k;
                #endif
                wi_flag = p->flag4(i,j,k);
                wi_cov = is_covered(p->level,i,j,k);   // always false now (covered cells skipped)
                if(projprobe)
                {
                    wi_div = V_lev*div; wi_Lp = Lp(i,j,k); wi_pc = pcorr(i,j,k);
                    wi_dax[0] = -V_lev*(u0(i,j,k) - u0(i-1,j,k))/(alpha*p->dt*p->DXN[IP]);
                    wi_dax[1] = p->j_dir ? -V_lev*(v0(i,j,k) - v0(i,j-1,k))/(alpha*p->dt*p->DYN[JP]) : 0.0;
                    wi_dax[2] = -V_lev*(w0(i,j,k) - w0(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                    wi_pcn[0]=pcorr(i-1,j,k); wi_pcn[1]=pcorr(i+1,j,k);
                    wi_pcn[2]=pcorr(i,j-1,k); wi_pcn[3]=pcorr(i,j+1,k);
                    wi_pcn[4]=pcorr(i,j,k-1); wi_pcn[5]=pcorr(i,j,k+1);
                    wi_fln[0]=p->flag4(i-1,j,k); wi_fln[1]=p->flag4(i+1,j,k);
                    wi_fln[2]=p->flag4(i,j-1,k); wi_fln[3]=p->flag4(i,j+1,k);
                    wi_fln[4]=p->flag4(i,j,k-1); wi_fln[5]=p->flag4(i,j,k+1);
                    wi_covn[0]=nb_covered(p->level,i-1,j,k); wi_covn[1]=nb_covered(p->level,i+1,j,k);
                    wi_covn[2]=nb_covered(p->level,i,j-1,k); wi_covn[3]=nb_covered(p->level,i,j+1,k);
                    wi_covn[4]=nb_covered(p->level,i,j,k-1); wi_covn[5]=nb_covered(p->level,i,j,k+1);
                    wi_rof[0]=pd->roface(p,a,-1,0,0); wi_rof[1]=pd->roface(p,a,1,0,0);
                    wi_rof[2]=pd->roface(p,a,0,-1,0); wi_rof[3]=pd->roface(p,a,0,1,0);
                    wi_rof[4]=pd->roface(p,a,0,0,-1); wi_rof[5]=pd->roface(p,a,0,0,1);
                    wi_phis=a->phi(i,j,k);
                    wi_phin[0]=a->phi(i-1,j,k); wi_phin[1]=a->phi(i+1,j,k);
                    wi_phin[2]=a->phi(i,j-1,k); wi_phin[3]=a->phi(i,j+1,k);
                    wi_phin[4]=a->phi(i,j,k-1); wi_phin[5]=a->phi(i,j,k+1);
                }
            }
        }
    }

    const double g_int = pgc->globalmax(r_int);
    const double g_bnd = pgc->globalmax(r_bnd);
    if(p->mpirank==0)
    std::cout<<"\n  [projcheck] max|L*pcorr + R(dU)|  interior="<<g_int
             <<"  boundary="<<g_bnd<<std::endl;
    if(wi_lev>=0 && wi_res==g_int)
    {
    std::cout<<"  [projcheck] worst interior res="<<wi_res<<" at lev="<<wi_lev
             <<" ("<<wi[0]<<","<<wi[1]<<","<<wi[2]<<")  flag4="<<wi_flag
             <<"  covered="<<wi_cov<<"  cf="<<cf_tag(wi_lev,wi[0],wi[1],wi[2])<<std::endl;
    if(projprobe)
    {
        std::cout<<"  [projprobe] Lp="<<wi_Lp<<"  V*div="<<wi_div<<"  (res=|Lp+V*div|)"
                 <<"  pcorr c="<<wi_pc<<std::endl;
        std::cout<<"  [projprobe] V*div split  x="<<wi_dax[0]<<"  y="<<wi_dax[1]
                 <<"  z="<<wi_dax[2]<<"  (residual axis = the covered/C-F face)"<<std::endl;
        std::cout<<"  [projprobe] pcorr  x-/x+="<<wi_pcn[0]<<"/"<<wi_pcn[1]
                 <<"  y-/y+="<<wi_pcn[2]<<"/"<<wi_pcn[3]
                 <<"  z-/z+="<<wi_pcn[4]<<"/"<<wi_pcn[5]<<std::endl;
        std::cout<<"  [projprobe] flag4  x-/x+="<<wi_fln[0]<<"/"<<wi_fln[1]
                 <<"  y-/y+="<<wi_fln[2]<<"/"<<wi_fln[3]
                 <<"  z-/z+="<<wi_fln[4]<<"/"<<wi_fln[5]<<std::endl;
        std::cout<<"  [projprobe] covered x-/x+="<<wi_covn[0]<<"/"<<wi_covn[1]
                 <<"  y-/y+="<<wi_covn[2]<<"/"<<wi_covn[3]
                 <<"  z-/z+="<<wi_covn[4]<<"/"<<wi_covn[5]
                 <<"  (a covered neighbour => patch-adjacent face)"<<std::endl;
        // High precision: cout carries a sticky setprecision(3) from the piter/umax prints, which
        // rounds roface (~998, but in-band it departs from W1 by O(0.1-10)) to a flat "998" and
        // hides the real face-to-face variation. Print roface AND W1-roface (the departure from
        // pure water) at full precision, then restore the stream precision.
        const std::streamsize oldprec = std::cout.precision();
        std::cout<<std::setprecision(10);
        std::cout<<"  [projprobe] roface x-/x+="<<wi_rof[0]<<"/"<<wi_rof[1]
                 <<"  y-/y+="<<wi_rof[2]<<"/"<<wi_rof[3]
                 <<"  z-/z+="<<wi_rof[4]<<"/"<<wi_rof[5]
                 <<"  (corrector face density; compare vs |L.coeff| on the C-F face)"<<std::endl;
        std::cout<<"  [projprobe] W1-roface x-/x+="<<(p->W1-wi_rof[0])<<"/"<<(p->W1-wi_rof[1])
                 <<"  y-/y+="<<(p->W1-wi_rof[2])<<"/"<<(p->W1-wi_rof[3])
                 <<"  z-/z+="<<(p->W1-wi_rof[4])<<"/"<<(p->W1-wi_rof[5])
                 <<"  (0 = pure water; nonzero = in the density band)"<<std::endl;
        std::cout<<"  [projprobe] phi c="<<wi_phis<<"  x-/x+="<<wi_phin[0]<<"/"<<wi_phin[1]
                 <<"  y-/y+="<<wi_phin[2]<<"/"<<wi_phin[3]
                 <<"  z-/z+="<<wi_phin[4]<<"/"<<wi_phin[5]<<std::endl;
        std::cout<<std::setprecision(oldprec);
    }
    }
    if(wb_lev>=0 && wb_res==g_bnd)
    std::cout<<"  [projcheck] worst boundary res="<<wb_res<<" at lev="<<wb_lev
             <<" ("<<wb[0]<<","<<wb[1]<<","<<wb[2]<<")  flag4="<<wb_flag
             <<"  cf="<<cf_tag(wb_lev,wb[0],wb[1],wb[2])<<std::endl;

    if(adjoint)
    {
        const double g_cf   = pgc->globalmax(r_cf);
        const double g_surf = pgc->globalmax(r_surf);
        if(p->mpirank==0)
        std::cout<<"  [adjoint]  max|L*pcorr + R(dU)|  C-F="<<g_cf
                 <<"  surface="<<g_surf<<"   (both must be ~solver tol)"<<std::endl;
        if(wc_lev>=0 && r_cf==g_cf)
        {
            std::cout<<"  [adjoint]  worst C-F res="<<r_cf<<" at lev="<<wc_lev
                     <<" ("<<wc[0]<<","<<wc[1]<<","<<wc[2]<<")  flag4="<<wc_flag<<std::endl;
            std::cout<<"  [adjoint]  Lp="<<wc_Lp<<"  V*div="<<wc_div<<"  pcorr="<<wc_pc<<std::endl;
            std::cout<<"  [adjoint]  V*div split  x="<<wc_dax[0]<<"  y="<<wc_dax[1]
                     <<"  z="<<wc_dax[2]<<"   (residual axis = the C-F / refluxed face)"<<std::endl;
            std::cout<<"  [adjoint]  flag4   x-/x+="<<wc_fln[0]<<"/"<<wc_fln[1]
                     <<"  y-/y+="<<wc_fln[2]<<"/"<<wc_fln[3]
                     <<"  z-/z+="<<wc_fln[4]<<"/"<<wc_fln[5]<<std::endl;
            std::cout<<"  [adjoint]  covered x-/x+="<<wc_covn[0]<<"/"<<wc_covn[1]
                     <<"  y-/y+="<<wc_covn[2]<<"/"<<wc_covn[3]
                     <<"  z-/z+="<<wc_covn[4]<<"/"<<wc_covn[5]
                     <<"   (covered neighbour => that face is the C-F/reflux face)"<<std::endl;
        }
        if(ws_lev>=0 && r_surf==g_surf)
        std::cout<<"  [adjoint]  worst surface res="<<r_surf<<" at lev="<<ws_lev
                 <<" ("<<ws[0]<<","<<ws[1]<<","<<ws[2]<<")  flag4="<<ws_flag<<std::endl;
    }
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

// REEF_CF_PROJECTION_GROUP member (6b) — u/v/wpgrad build the predictor pressure force
// (press(i+1)-press(i))/(DXP*roface). At the C-F interface press's ghost ring is filled by
// FillCoarseFineCellGhost (gcv 41), so the coarse and fine predictor gradients see the
// d_cf-consistent value; without it the fine block gets a spurious hydrostatic force. Keep
// DXP/roface consistent with the group (canonical note in hypre_ssamg::amr_cf_coefficients).
void pjm_corr::upgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    double dp = 0.0;
    const bool relPressure = p->Y9;
    const bool skip_covered = (std::getenv("REEF_SKIP_COVERED_PGRAD") != nullptr);

    // REEF_HYDRO_PROBE (horizontal x): at rest the horizontal predictor force must be ~0 EVERYWHERE
    // (no horizontal gravity, no horizontal hydrostatic pressure gradient). Any nonzero |F_imb| is a
    // pure artifact; its location (surface / C-F / deep) localises the horizontal injection driving umax.
    const bool hydro_probe = (std::getenv("REEF_HYDRO_PROBE") != nullptr) && relPressure;
    double hp_w=0.0; int hp_l=-1, hp_i[3]={-1,-1,-1}; double hp_rof=0,hp_pg=0,hp_gg=0,hp_ps=0,hp_pn=0;
    double hp_prs=0,hp_prn=0; int hp_cs=0,hp_cn=0;   // press(self,i+1) + covered flags
    double hp_deep=0.0; int hp_dl=-1, hp_di[3]={-1,-1,-1};

    ULOOP
    {
        #if USE_AMREX
        // skip the predictor pgrad on any face touching a covered cell (self OR neighbour): those
        // faces are C-F/fine-authoritative and overwritten by reflux, and the covered press is
        // horizontally inconsistent (per-column hydrostatic offset). Neighbour term added 2026-07-06.
        if(skip_covered && (cell_is_covered(p,i,j,k) || cell_is_covered(p,i+1,j,k))) continue;
        #endif
        dp = a->press(i+1,j,k)-a->press(i,j,k);
        // if(relPressure) dp -= a->press0(i+1,j,k)-a->press0(i,j,k);
        a->F(i,j,k) -= PORVAL1*dp/(p->DXP[IP]*pd->roface(p,a,1,0,0));
        if(relPressure) a->F(i,j,k) += PORVAL1*(a->grav_pot(i+1,j,k)-a->grav_pot(i,j,k))/p->DXP[IP];

        if(hydro_probe)
        {
            #if USE_AMREX
            const double rof = pd->roface(p,a,1,0,0);
            const double pg  = dp/(p->DXP[IP]*rof);
            const double gg  = (a->grav_pot(i+1,j,k)-a->grav_pot(i,j,k))/p->DXP[IP];
            const double imb = std::fabs(PORVAL1*(-pg + gg));
            const double phis=a->phi(i,j,k), phin=a->phi(i+1,j,k);
            if(imb > hp_w)
            {
                hp_w=imb; hp_l=p->level;
                hp_i[0]=i+p->amr_tile_lo.x; hp_i[1]=j+p->amr_tile_lo.y; hp_i[2]=k+p->amr_tile_lo.z;
                hp_rof=rof; hp_pg=pg; hp_gg=gg; hp_ps=phis; hp_pn=phin;
                hp_prs=a->press(i,j,k); hp_prn=a->press(i+1,j,k);
                #if USE_AMREX
                hp_cs=cell_is_covered(p,i,j,k)?1:0; hp_cn=cell_is_covered(p,i+1,j,k)?1:0;
                #endif
            }
            if(std::min(std::fabs(phis),std::fabs(phin)) > 0.1 && imb > hp_deep)
            {
                hp_deep=imb; hp_dl=p->level;
                hp_di[0]=i+p->amr_tile_lo.x; hp_di[1]=j+p->amr_tile_lo.y; hp_di[2]=k+p->amr_tile_lo.z;
            }
            #endif
        }
    }

    if(hydro_probe && p->mpirank==0)
    {
        std::cout<<"  [hydroprobe-u] worst |F_imb|="<<hp_w<<" at lev="<<hp_l
                 <<" ("<<hp_i[0]<<","<<hp_i[1]<<","<<hp_i[2]<<")  roface="<<hp_rof
                 <<" pgrad/rof="<<hp_pg<<" gpot_grad="<<hp_gg
                 <<" phi[self,i+1]="<<hp_ps<<","<<hp_pn
                 <<" press[self,i+1]="<<hp_prs<<","<<hp_prn
                 <<" covered[self,i+1]="<<hp_cs<<","<<hp_cn<<std::endl;
        std::cout<<"  [hydroprobe-u] worst DEEP |F_imb|="<<hp_deep<<" at lev="<<hp_dl
                 <<" ("<<hp_di[0]<<","<<hp_di[1]<<","<<hp_di[2]<<")  (should be ~0)"<<std::endl;
    }
}

void pjm_corr::vpgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    if(p->j_dir)
    {
        double dp = 0.0;
        const bool relPressure = p->Y9;
        const bool skip_covered = (std::getenv("REEF_SKIP_COVERED_PGRAD") != nullptr);

        // REEF_HYDRO_PROBE (horizontal y): mirror of upgrad; at rest |G_imb| must be ~0 everywhere.
        const bool hydro_probe = (std::getenv("REEF_HYDRO_PROBE") != nullptr) && relPressure;
        double hp_w=0.0; int hp_l=-1, hp_i[3]={-1,-1,-1}; double hp_rof=0,hp_pg=0,hp_gg=0,hp_ps=0,hp_pn=0;
        double hp_deep=0.0; int hp_dl=-1, hp_di[3]={-1,-1,-1};

        VLOOP
        {
            #if USE_AMREX
            if(skip_covered && (cell_is_covered(p,i,j,k) || cell_is_covered(p,i,j+1,k))) continue;
            #endif
            dp = a->press(i,j+1,k)-a->press(i,j,k);
            // if(relPressure) dp -= a->press0(i,j+1,k)-a->press0(i,j,k);
            a->G(i,j,k) -= PORVAL2*dp/(p->DYP[JP]*pd->roface(p,a,0,1,0));
            if(relPressure) a->G(i,j,k) += PORVAL2*(a->grav_pot(i,j+1,k)-a->grav_pot(i,j,k))/p->DYP[JP];

            if(hydro_probe)
            {
                #if USE_AMREX
                const double rof = pd->roface(p,a,0,1,0);
                const double pg  = dp/(p->DYP[JP]*rof);
                const double gg  = (a->grav_pot(i,j+1,k)-a->grav_pot(i,j,k))/p->DYP[JP];
                const double imb = std::fabs(PORVAL2*(-pg + gg));
                const double phis=a->phi(i,j,k), phin=a->phi(i,j+1,k);
                if(imb > hp_w)
                {
                    hp_w=imb; hp_l=p->level;
                    hp_i[0]=i+p->amr_tile_lo.x; hp_i[1]=j+p->amr_tile_lo.y; hp_i[2]=k+p->amr_tile_lo.z;
                    hp_rof=rof; hp_pg=pg; hp_gg=gg; hp_ps=phis; hp_pn=phin;
                }
                if(std::min(std::fabs(phis),std::fabs(phin)) > 0.1 && imb > hp_deep)
                {
                    hp_deep=imb; hp_dl=p->level;
                    hp_di[0]=i+p->amr_tile_lo.x; hp_di[1]=j+p->amr_tile_lo.y; hp_di[2]=k+p->amr_tile_lo.z;
                }
                #endif
            }
        }

        if(hydro_probe && p->mpirank==0)
        {
            std::cout<<"  [hydroprobe-v] worst |G_imb|="<<hp_w<<" at lev="<<hp_l
                     <<" ("<<hp_i[0]<<","<<hp_i[1]<<","<<hp_i[2]<<")  roface="<<hp_rof
                     <<" pgrad/rof="<<hp_pg<<" gpot_grad="<<hp_gg
                     <<" phi[self,j+1]="<<hp_ps<<","<<hp_pn<<std::endl;
            std::cout<<"  [hydroprobe-v] worst DEEP |G_imb|="<<hp_deep<<" at lev="<<hp_dl
                     <<" ("<<hp_di[0]<<","<<hp_di[1]<<","<<hp_di[2]<<")  (should be ~0)"<<std::endl;
        }
    }
}

void pjm_corr::wpgrad(lexer*p, fdm* a, slice &eta, slice &eta_n)
{
    double dp = 0.0;
    const bool relPressure = p->Y9;

    // REEF_HYDRO_PROBE: at-rest vertical hydrostatic imbalance = wpgrad's contribution to H.
    // = -grad_z(press)/roface + grad_z(grav_pot); must be ~0 everywhere when well-balanced.
    // Report the global-worst cell + its components/phi, AND the worst among DEEP cells (away
    // from the surface). If the imbalance concentrates at surface cells on the fine patch and is
    // ~0 deep, the well-balancing breaks specifically at the surface-crossing C-F interface.
    const bool hydro_probe = (std::getenv("REEF_HYDRO_PROBE") != nullptr) && relPressure;
    const bool skip_covered = (std::getenv("REEF_SKIP_COVERED_PGRAD") != nullptr);
    double hp_w=0.0; int hp_l=-1, hp_i[3]={-1,-1,-1}; double hp_rof=0,hp_pg=0,hp_gg=0,hp_ps=0,hp_pn=0;
    double hp_deep=0.0; int hp_dl=-1, hp_di[3]={-1,-1,-1};

    WLOOP
    {
        #if USE_AMREX
        if(skip_covered && (cell_is_covered(p,i,j,k) || cell_is_covered(p,i,j,k+1))) continue;
        #endif
        dp = a->press(i,j,k+1)-a->press(i,j,k);
        // if(relPressure) dp -= a->press0(i,j,k+1)-a->press0(i,j,k);
        a->H(i,j,k) -= PORVAL3*dp/(p->DZP[KP]*pd->roface(p,a,0,0,1));
        if(relPressure) a->H(i,j,k) += PORVAL3*(a->grav_pot(i,j,k+1)-a->grav_pot(i,j,k))/p->DZP[KP];

        if(hydro_probe)
        {
            #if USE_AMREX
            const double rof = pd->roface(p,a,0,0,1);
            const double pg  = dp/(p->DZP[KP]*rof);
            const double gg  = (a->grav_pot(i,j,k+1)-a->grav_pot(i,j,k))/p->DZP[KP];
            const double imb = std::fabs(PORVAL3*(-pg + gg));
            const double phis=a->phi(i,j,k), phin=a->phi(i,j,k+1);
            if(imb > hp_w)
            {
                hp_w=imb; hp_l=p->level;
                hp_i[0]=i+p->amr_tile_lo.x; hp_i[1]=j+p->amr_tile_lo.y; hp_i[2]=k+p->amr_tile_lo.z;
                hp_rof=rof; hp_pg=pg; hp_gg=gg; hp_ps=phis; hp_pn=phin;
            }
            // DEEP = both faces well away from the surface (|phi| > 0.1 m ~ 2 coarse cells)
            if(std::min(std::fabs(phis),std::fabs(phin)) > 0.1 && imb > hp_deep)
            {
                hp_deep=imb; hp_dl=p->level;
                hp_di[0]=i+p->amr_tile_lo.x; hp_di[1]=j+p->amr_tile_lo.y; hp_di[2]=k+p->amr_tile_lo.z;
            }
            #endif
        }
    }

    if(hydro_probe && p->mpirank==0)
    {
        std::cout<<"  [hydroprobe-w] worst |H_imb|="<<hp_w<<" at lev="<<hp_l
                 <<" ("<<hp_i[0]<<","<<hp_i[1]<<","<<hp_i[2]<<")  roface="<<hp_rof
                 <<" pgrad/rof="<<hp_pg<<" gpot_grad="<<hp_gg
                 <<" phi[self,k+1]="<<hp_ps<<","<<hp_pn<<std::endl;
        std::cout<<"  [hydroprobe-w] worst DEEP |H_imb|="<<hp_deep<<" at lev="<<hp_dl
                 <<" ("<<hp_di[0]<<","<<hp_di[1]<<","<<hp_di[2]<<")  (should be ~0)"<<std::endl;
    }
}

void pjm_corr::ini(lexer*p, fdm* a, ghostcell *pgc)
{
    reference_ini(p,a,pgc);
}

// Composite predictor divergence (D side): make the PREDICTOR normal velocity single-valued
// at the C-F interface -- set each covered coarse face = area-average of the fine faces across
// it, so the coarse cell's divergence (rhs) uses the summed fine flux, exactly the flux L's
// coarse row couples to (conservative M_cf*V_c = M_fc*V_f). This makes the discrete divergence D
// the adjoint of the gradient G the matrix used, so D(1/rho)G = L. No-op for nlevs<=1.
void pjm_corr::cf_predictor_sync(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel, solver* psolv, double alpha)
{
    #if USE_AMREX
    if(p->nlevs > 1)
    {
        // C-F-aware predictor: slave the fine HIGH C-F normal faces to the clean coarse predictor
        // velocity BEFORE the reflux, so the predictor never hands the projection a leaked C-F
        // velocity it cannot fully remove. psolv here is ppoissonsolv (hypre_ssamg), which owns
        // cf_links. Gated on Y9 so it toggles with the well-balanced setup under test.
        if(p->Y9)
        psolv->cf_velocity_fill_from_coarse(p,a,pgc,uvel,vvel,wvel);

        cf_average_down_velocity(p,uvel,vvel,wvel);

        // average_down updated COARSE valid cells under the fine patch; refresh the halos so a
        // neighbour rank's rhs reads the single-valued interface velocity (not the stale
        // pre-average value).
        vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);
    }
    #endif
}

// B5 (PLAN_Rebalance_PhiSync.md): re-equilibrate press against the CURRENT density after reinit
// moved phi. The predictor buoyancy force b = -grad(press)/rho_f + grad(grav_pot) (Y9 form, the
// same expressions as u/v/wpgrad) is no longer zero because press still matches the pre-reinit
// density. b splits into a curl-free part (removable by a pressure re-solve) and a rotational
// part (grad(1/rho) x grad(p), unreachable by any pressure). One projection of the pure-buoyancy
// predictor from rest -- utmp = alpha*dt*b -- and keeping only the pressure update (press +=
// pcorr, NO velcorr) removes the curl-free part to solver tolerance. The rotational remainder is
// the A-family follow-up. Velocity is untouched: this is a pressure re-initialisation, not a
// time step. Structurally identical to start() with the "velocity" being the buoyancy source.
void pjm_corr::rebalance(lexer*p, fdm* a, ghostcell *pgc, poisson* ppois, solver* psolv, ioflow* pflow)
{
    if(std::getenv("REEF_REBALANCE") == nullptr)
    return;

    // Y9 (relative-pressure / well-balanced) form only in this pass. The non-Y9 form is the same
    // construction with the body-force gravity a->gi/gj/gk (added in momentum irhs/jrhs/krhs)
    // replacing the grav_pot gradient. TODO: implement the non-Y9 branch.
    if(!p->Y9)
    {
        if(p->mpirank==0 && p->count<=1)
        cout<<"  [rebalance] skipped: only the Y9 (relPressure) form is implemented"<<endl;
        return;
    }

    const double starttime = pgc->timer();

    // Optional extra passes for measurement only (REEF_REBALANCE_PASSES, default 1). Pass 1
    // removes the entire curl-free part; further passes only expose the rotational floor.
    const int passes = (std::getenv("REEF_REBALANCE_PASSES") ? std::atoi(std::getenv("REEF_REBALANCE_PASSES")) : 1);

    // alpha cancels: utmp = alpha*dt*b and rhs() divides the divergence by alpha*dt, so the RHS
    // is div(b) independent of alpha. Use 1.0 -- rebalance is a standalone re-init, not an RK stage.
    const double alpha = 1.0;

    // Face temporaries for the buoyancy predictor. Self-allocating (own MultiFab), like the
    // projection-consistency probe's field4 locals; do NOT reuse a->F/G/H (momentum-stage state).
    field1 utmp(p);
    field2 vtmp(p);
    field3 wtmp(p);

    for(int pass=0; pass<passes; ++pass)
    {
        utmp.setVal(0.0,true);
        vtmp.setVal(0.0,true);
        wtmp.setVal(0.0,true);

        // b at rest = the u/v/wpgrad force (grad(press)/rho_f + grav_pot gradient). Kept in lockstep
        // with those methods; utmp = alpha*dt*b is exactly the at-rest predictor with no convec/diff.
        ULOOP
        utmp(i,j,k) = alpha*p->dt*PORVAL1*
            ( -(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*pd->roface(p,a,1,0,0))
              +(a->grav_pot(i+1,j,k)-a->grav_pot(i,j,k))/p->DXP[IP] );

        if(p->j_dir==1)
        VLOOP
        vtmp(i,j,k) = alpha*p->dt*PORVAL2*
            ( -(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*pd->roface(p,a,0,1,0))
              +(a->grav_pot(i,j+1,k)-a->grav_pot(i,j,k))/p->DYP[JP] );

        WLOOP
        wtmp(i,j,k) = alpha*p->dt*PORVAL3*
            ( -(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1))
              +(a->grav_pot(i,j,k+1)-a->grav_pot(i,j,k))/p->DZP[KP] );

        // Same pipeline as start(), minus velcorr: halo fill -> composite C-F sync -> divergence
        // rhs -> assemble matrix with the CURRENT roface -> solve -> press += pcorr.
        vel_setup(p,a,pgc,utmp,vtmp,wtmp,alpha);
        cf_predictor_sync(p,a,pgc,utmp,vtmp,wtmp,psolv,alpha);
        rhs(p,a,pgc,utmp,vtmp,wtmp,alpha);   // also zeroes pcorr

        ppois->start(p,a,pcorr);
        psolv->start(p,a,pgc,pcorr,a->rhsvec,5);

        #if USE_AMREX
        const int gcval_press = (p->nlevs > 1) ? 41 : 40;
        #else
        const int gcval_press = 40;
        #endif
        pgc->start4(p,pcorr,gcval_press,false);

        presscorr(p,a);                      // press += pcorr  (NO velcorr)
        reference_start(p,a,pgc);
        pgc->start4(p,a->press,gcval_press,false);
    }

    // Not added to p->poissontime -- momentum's start() overwrites it this step; report separately.
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"  [rebalance] piter: "<<p->solveriter<<"  pres: "<<setprecision(3)<<p->final_res
        <<"  ptime: "<<setprecision(3)<<pgc->timer()-starttime<<endl;
}
