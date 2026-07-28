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

#include"levelset_RK3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
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
#include"heaviside_ls.h"
#include<cstdlib>
#include<cmath>
#if USE_AMREX
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#endif

levelset_RK3::levelset_RK3(lexer* p, fdm *a, ghostcell* pgc, heat *&pheat, concentration *&pconc):gradient(p), ark1(p), ark2(p)
{
    if(p->F50==1)
    gcval_phi=51;
    else if(p->F50==2)
    gcval_phi=52;
    else if(p->F50==3)
    gcval_phi=53;
    else if(p->F50==4)
    gcval_phi=54;

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
    if(p->F300>0)
    pupdate = new fluid_update_void();

    if(p->F46==2)
    ppicard = new picard_f(p);
    else if(p->F46==3)
    ppicard = new picard_lsm(p);
    else
    ppicard = new picard_void(p);

    gcval_u=10;
    gcval_v=11;
    gcval_w=12;
}

levelset_RK3::~levelset_RK3()
{
    delete pupdate;
    delete ppicard;
}

void levelset_RK3::start(fdm* a, lexer* p, convection* pconvec, solver*, ghostcell* pgc, ioflow* pflow, reini* preini, field &ls)
{
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkin(p,a,pgc,ark2);
    pflow->fsfrkout(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark2);
    ppicard->volcalc(p,a,pgc,ls);

    const double dt = p->dt;

    double setValtime, pconvectime, calctime, phirelaxtime, gctime, solidforcetime, reinitime, picardtime, update_time;
    double temptime;

// Step 1
    starttime=pgc->timer();

    a->L.setVal(0.0);

    setValtime=pgc->timer();

    pconvec->start(p,a,ls,4,a->u,a->v,a->w);

    pconvectime=pgc->timer();
    FIELDLOOP(
        ark1,
        FIELD_CONST(ls); FIELD_CONST_MEMBER(a, L),
        ark1(i,j,k) = ls(i,j,k) + dt*member_L(i,j,k);
    )
    calctime=pgc->timer();

    pflow->phi_relax(p,pgc,ark1);

    phirelaxtime=pgc->timer();

    pgc->start4(p,ark1,gcval_phi);
    gctime=pgc->timer();

    pgc->solid_forcing_lsm(p,a,ark1);
    solidforcetime=pgc->timer();

// Step 2
    a->L.setVal(0.0);

    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);

    FIELDLOOP(
        ark2,
        FIELD_CONST(ls); FIELD_CONST(ark1); FIELD_CONST_MEMBER(a, L),
        ark2(i,j,k) = 0.75*ls(i,j,k) + 0.25*ark1(i,j,k) + 0.25*dt*member_L(i,j,k);
    )

    pflow->phi_relax(p,pgc,ark2);

    pgc->start4(p,ark2,gcval_phi);
    pgc->solid_forcing_lsm(p,a,ark2);

// Step 3
    a->L.setVal(0.0);

    pconvec->start(p,a,ark2,4,a->u,a->v,a->w);

    FIELDLOOP(
        ls,
        FIELD_CONST(ark2) FIELD_CONST_MEMBER(a, L),
        ls(i,j,k) = (1.0/3.0)*ls(i,j,k) + (2.0/3.0)*ark2(i,j,k) + (2.0/3.0)*dt*member_L(i,j,k);
    )

    pflow->phi_relax(p,pgc,ls);

    pgc->start4(p,ls,gcval_phi);
    pgc->solid_forcing_lsm(p,a,ls);

    p->lsmtime=pgc->timer()-starttime;

    temptime=pgc->timer();

    preini->start(a,p,ls,pgc,pflow);

    reinitime=pgc->timer();

    ppicard->correct_ls(p,a,pgc,ls);
    picardtime=pgc->timer();

    #if USE_AMREX
    if(p->nlevs > 1)
    {
        // No explicit phi average_down here: ghostcell::start4's AMReX path average_downs every
        // cell field on each call (gc_start.cpp, do_avgdown default), including reini's own
        // per-substep start4 calls — covered coarse phi is already the fine average. This start4
        // only refreshes halos/C-F ghosts after correct_ls/picardmove touched phi.
        pgc->start4(p, ls, gcval_phi);

        if(std::getenv("REEF_PHI_CF_PROBE"))
            phi_cf_probe(p, a, pgc, ls);
    }
    #endif

	pupdate->start(p,a,pgc,a->u,a->v,a->w);
    update_time=pgc->timer();
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"lsmtime: "<<setprecision(3)<<p->lsmtime<<endl;

    if(p->mpirank==0 && p->count%p->P12==0 && p->count>0 && std::getenv("REEF_timing"))
    {
        const int precision = 5;
        std::cout<<"\n"
        <<"Timing for level set RK3 step "<<p->count<<"\n"
        <<"\tsteptime:     "<<std::setprecision(precision)<<gctime-starttime<<"\n"
        <<"\tsetValtime:  "<<std::setprecision(precision)<<setValtime-starttime<<"\n"
        <<"\tpconvectime: "<<std::setprecision(precision)<<pconvectime-setValtime<<"\n"
        <<"\tcalctime:    "<<std::setprecision(precision)<<calctime-pconvectime<<"\n"
        <<"\tgctime:      "<<std::setprecision(precision)<<gctime-phirelaxtime<<"\n\n"
        <<"\treinitime:   "<<std::setprecision(precision)<<reinitime-temptime<<"\n"
        <<"\tupdate_time: "<<std::setprecision(precision)<<update_time-picardtime<<"\n"<<std::endl;
    }
}

void levelset_RK3::update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    pupdate->start(p,a,pgc,a->u,a->v,a->w);
}

#if USE_AMREX
// REEF_PHI_CF_PROBE (Plan Part 1 / C7) ----------------------------------------------------------
// Verify single-valued composite phi/density at the C-F interface AFTER the post-reinit sync
// (start4 average_downs phi internally). All loops are raw AMReX MFIter/LoopOnCpu in GLOBAL
// indices — no amr_tile_lo translation (that offset only applies to the LOOP-macro accessors).
// Three metrics, per adjacent level pair (clev, flev):
//
//   (A) max |avg(phi_fine) - phi_coarse| over covered coarse cells. Plumbing check for the
//       start4-internal average_down; the copy leaves uncovered cells self-equal (diff 0).
//       Expect ~machine eps.
//
//   (B) max |roface_fine - roface_coarse| over covered coarse faces, split by face class
//       (covered mask = amr_cell_mf: 1 uncovered, 0 covered):
//         ring DEEP    — exactly one adjacent cell covered (true C-F interface face) AND the
//                        coarse Heaviside saturated (|phi_mean| >= psi): both sides must be the
//                        identical constant W1 or W3, so ANY nonzero value is a hard plumbing
//                        bug (ghost fill / indexing / BC). This is the pass/fail line.
//         ring in-band — C-F face inside the Heaviside band: upper bound only. Contains the
//                        irreducible two-level part (transverse phi variation across the r^2
//                        sub-faces + Heaviside nonlinearity), NOT attributable to the ghost
//                        fill; metric (C) isolates the fill itself.
//         interior     — both cells covered: coarse side is H(avg phi), fine side avg(H(phi)),
//                        differing by the Jensen gap of the nonlinear H. Fine-authoritative
//                        region, solver-irrelevant; reported separately so it cannot mask the
//                        ring signal.
//       Domain-boundary faces are skipped so the ring metric is pure C-F. roface replicates
//       density_f::roface exactly (H(0.5*(phi_a+phi_b), p->psi); psi is a global constant).
//
//   (C) fine C-F ghost phi: current fill (FillPatchTwoLevels / cell_cons_interp) vs the matrix
//       d_cf convention (FillCoarseFineCellGhost pass 1: g = b + 2/(1+r)*(c - b)). Reported as
//       max|dphi| and — decisively — max|droface|, the change in the fine-face density the
//       predictor would see if phi got the gcv-41-style fill. ~0 => the cell_cons fill is
//       density-adequate and C7 closes; O(W1-W3) => extend the matrix-consistent fill to the
//       phi gcvals (plan Part 1, item 2).
void levelset_RK3::phi_cf_probe(lexer* p, fdm* a, ghostcell* pgc, field& ls)
{
    if(p->nlevs <= 1) return;

    const double psi = p->psi;
    const double W1 = p->W1, W3 = p->W3;

    auto roface = [&] (double phia, double phib) -> double
    {
        const double H = heaviside_ls(0.5*(phia+phib), psi);
        return W1*H + W3*(1.0-H);
    };

    double phi_max  = 0.0;                                    // (A)
    double rof_deep = 0.0, rof_band = 0.0, rof_int = 0.0;     // (B)
    int    deep_loc[5] = {-1,0,0,0,0}, band_loc[5] = {-1,0,0,0,0};
    double gphi_max = 0.0, grof_max = 0.0;                    // (C)
    int    gc_loc[5] = {-1,0,0,0,0};

    for(int clev = 0; clev <= p->nlevs-2; ++clev)
    {
        const int flev = clev + 1;

        // (A) covered coarse phi vs fine average --------------------------------------------
        amrex::MultiFab phi_avg(p->amrex_box_array[clev],
                                p->amrex_distribution_mapping[clev], 1, 0);
        amrex::MultiFab::Copy(phi_avg, ls.GetMultiFab(clev), 0, 0, 1, 0);
        amrex::average_down(ls.GetMultiFab(flev), phi_avg, 0, 1, p->ref_vec);

        for(amrex::MFIter mfi(phi_avg); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const auto pa = phi_avg.const_array(mfi);
            const auto pc = ls.GetMultiFab(clev).const_array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
            { phi_max = std::max(phi_max, std::fabs(pa(i,j,k) - pc(i,j,k))); });
        }

        // Covered mask with a 1-cell ghost so the face-adjacent cell's state is known across
        // box boundaries (same construction as FillCoarseFineCellGhost pass 2). Domain-exterior
        // ghosts stay 1 (uncovered) — combined with the cdom check below, boundary faces drop out.
        amrex::iMultiFab maskg(p->amr_cell_mf[clev].boxArray(),
                               p->amr_cell_mf[clev].DistributionMap(), 1, 1);
        maskg.setVal(1);
        amrex::iMultiFab::Copy(maskg, p->amr_cell_mf[clev], 0, 0, 1, 0);
        maskg.FillBoundary(p->amrex_geometry[clev].periodicity());
        const amrex::Box cdom = p->amrex_geometry[clev].Domain();

        // (B) covered C-F face roface: fine-side average vs coarse-side ----------------------
        for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if(dir == 1 && p->j_dir != 1) continue;   // degenerate y in a 2D x-z run

            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

            amrex::MultiFab rf_fine(amrex::convert(p->amrex_box_array[flev], e),
                                    p->amrex_distribution_mapping[flev], 1, 0);
            amrex::MultiFab rf_crse(amrex::convert(p->amrex_box_array[clev], e),
                                    p->amrex_distribution_mapping[clev], 1, 0);

            // roface on the faces of each level (face f sits between cells f-e and f)
            auto fill_rof = [&] (amrex::MultiFab& rmf, int amrlev)
            {
                const auto& phimf = ls.GetMultiFab(amrlev);
                for(amrex::MFIter mfi(rmf); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& fbx = mfi.validbox();
                    auto       r  = rmf.array(mfi);
                    const auto ph = phimf.const_array(mfi);
                    amrex::LoopOnCpu(fbx, [&] (int i, int j, int k) noexcept
                    { r(i,j,k) = roface(ph(i-e[0], j-e[1], k-e[2]), ph(i,j,k)); });
                }
            };
            fill_rof(rf_fine, flev);
            fill_rof(rf_crse, clev);

            // copy coarse roface, then overwrite covered faces with the fine-face average;
            // uncovered faces stay self-equal so the diff isolates the covered faces.
            amrex::MultiFab rf_from_fine(rf_crse.boxArray(), rf_crse.DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(rf_from_fine, rf_crse, 0, 0, 1, 0);
            amrex::average_down_faces(rf_fine, rf_from_fine, p->ref_vec, p->amrex_geometry[clev]);

            for(amrex::MFIter mfi(rf_crse); mfi.isValid(); ++mfi)
            {
                const amrex::Box& fbx = mfi.validbox();
                const auto rc = rf_crse.const_array(mfi);
                const auto rf = rf_from_fine.const_array(mfi);
                const auto ph = ls.GetMultiFab(clev).const_array(mfi);
                const auto mk = maskg.const_array(mfi);
                amrex::LoopOnCpu(fbx, [&] (int i, int j, int k) noexcept
                {
                    const double d = std::fabs(rf(i,j,k) - rc(i,j,k));
                    if(d == 0.0) return;   // uncovered face (untouched copy) or exact match

                    const amrex::IntVect chi(i,j,k);
                    const amrex::IntVect clo = chi - e;
                    if(!cdom.contains(clo) || !cdom.contains(chi)) return;   // domain-boundary face

                    const bool cov_lo = (mk(clo) == 0);
                    const bool cov_hi = (mk(chi) == 0);

                    if(cov_lo && cov_hi)               // interior covered: Jensen-gap noise
                    { rof_int = std::max(rof_int, d); return; }
                    if(cov_lo == cov_hi) return;       // both uncovered: not written (defensive)

                    if(std::fabs(0.5*(ph(clo) + ph(chi))) >= psi)   // coarse H saturated
                    {
                        if(d > rof_deep)
                        { rof_deep = d; deep_loc[0]=clev; deep_loc[1]=dir; deep_loc[2]=i; deep_loc[3]=j; deep_loc[4]=k; }
                    }
                    else if(d > rof_band)
                    { rof_band = d; band_loc[0]=clev; band_loc[1]=dir; band_loc[2]=i; band_loc[3]=j; band_loc[4]=k; }
                });
            }
        }

        // (C) fine C-F ghost ring: current fill vs the d_cf convention ------------------------
        {
            const amrex::MultiFab& fmf  = ls.GetMultiFab(flev);
            const amrex::Box       fdom = p->amrex_geometry[flev].Domain();
            const amrex::BoxArray& fba  = fmf.boxArray();

            // Coarse phi on the coarsened-fine layout with a 2-cell ghost so coarsen(ghost) is
            // always addressable (mirrors FillCoarseFineCellGhost pass 1).
            amrex::MultiFab cphi(amrex::coarsen(fba, p->ref_vec), fmf.DistributionMap(), 1, 2);
            cphi.setVal(0.0);
            cphi.ParallelCopy(ls.GetMultiFab(clev), 0, 0, 1,
                              amrex::IntVect(0), amrex::IntVect(2),
                              p->amrex_geometry[clev].periodicity());

            for(amrex::MFIter mfi(fmf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& vbx = mfi.validbox();
                const amrex::Box& fbx = mfi.fabbox();
                const auto fa = fmf.const_array(mfi);
                const auto ca = cphi.const_array(mfi);

                for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                {
                    if(dir == 1 && p->j_dir != 1) continue;

                    const double         s = 2.0/(1.0 + double(p->ref_vec[dir]));
                    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

                    for(int side = 0; side < 2; ++side)
                    {
                        // first ghost layer outboard of the valid box in `dir`
                        amrex::Box gbx = vbx;
                        if(side == 0) { gbx.setSmall(dir, vbx.smallEnd(dir)-1); gbx.setBig(dir, vbx.smallEnd(dir)-1); }
                        else          { gbx.setSmall(dir, vbx.bigEnd(dir)+1);   gbx.setBig(dir, vbx.bigEnd(dir)+1);   }
                        gbx &= fbx;
                        if(gbx.isEmpty()) continue;

                        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) noexcept
                        {
                            const amrex::IntVect g(i,j,k);
                            if(!fdom.contains(g)) return;   // domain-boundary ghost: BC owns it
                            if(fba.contains(g))   return;   // fine-fine ghost: FillBoundary owns it

                            const amrex::IntVect b = (side == 0) ? g + e : g - e;   // interior boundary cell
                            const double phb    = fa(b);
                            const double phc    = ca(amrex::coarsen(g, p->ref_vec));
                            const double ph_dcf = phb + s*(phc - phb);
                            const double ph_cur = fa(i,j,k);

                            gphi_max = std::max(gphi_max, std::fabs(ph_cur - ph_dcf));

                            const double dr = std::fabs(roface(ph_cur, phb) - roface(ph_dcf, phb));
                            if(dr > grof_max)
                            { grof_max = dr; gc_loc[0]=flev; gc_loc[1]=dir; gc_loc[2]=i; gc_loc[3]=j; gc_loc[4]=k; }
                        });
                    }
                }
            }
        }
    }

    double glob[6] = {phi_max, rof_deep, rof_band, rof_int, gphi_max, grof_max};
    amrex::ParallelDescriptor::ReduceRealMax(glob, 6);

    if(p->mpirank==0)
    {
        std::cout<<"  [phi_cf_probe] (A) covered-cell max|avg(phi_f)-phi_c| = "<<glob[0]
                 <<"   (plumbing check; expect ~1e-16)"<<std::endl;
        std::cout<<"  [phi_cf_probe] (B) roface ring DEEP        = "<<glob[1]
                 <<"   (MUST be 0; nonzero = ghost/indexing/BC bug)"<<std::endl;
        std::cout<<"  [phi_cf_probe] (B) roface ring in-band     = "<<glob[2]
                 <<"   (upper bound; contains irreducible two-level part)"<<std::endl;
        std::cout<<"  [phi_cf_probe] (B) roface interior covered = "<<glob[3]
                 <<"   (Jensen-gap noise; solver-irrelevant)"<<std::endl;
        std::cout<<"  [phi_cf_probe] (C) ghost fill vs d_cf: max|dphi| = "<<glob[4]
                 <<"  max|droface| = "<<glob[5]
                 <<"   (decides gcv-41-on-phi; ~0 closes C7)"<<std::endl;
    }

    // worst locations, printed by the owning rank (same pattern as projcheck)
    if(rof_deep > 0.0 && rof_deep == glob[1])
    std::cout<<"  [phi_cf_probe] worst ring DEEP at clev="<<deep_loc[0]<<" dir="<<deep_loc[1]
             <<" face("<<deep_loc[2]<<","<<deep_loc[3]<<","<<deep_loc[4]<<")"<<std::endl;
    if(rof_band > 0.0 && rof_band == glob[2])
    std::cout<<"  [phi_cf_probe] worst ring in-band at clev="<<band_loc[0]<<" dir="<<band_loc[1]
             <<" face("<<band_loc[2]<<","<<band_loc[3]<<","<<band_loc[4]<<")"<<std::endl;
    if(grof_max > 0.0 && grof_max == glob[5])
    std::cout<<"  [phi_cf_probe] worst (C) droface at flev="<<gc_loc[0]<<" dir="<<gc_loc[1]
             <<" ghost("<<gc_loc[2]<<","<<gc_loc[3]<<","<<gc_loc[4]<<")"<<std::endl;
}
#endif
