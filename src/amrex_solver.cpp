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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#include "amrex_solver.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "field1.h"
#include "field2.h"
#include "field3.h"
#include "field4.h"
#include "density_f.h"
#include "density_pst.h"
#include "density_comp.h"
#include "density_vof.h"
#include "density_rheo.h"
#include "looping.h"
#include "definitions_amrex.h"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>

#include <iostream>
#include <iomanip>
#include <string>

// REEF_MLMG_FPE: scoped invalid-FP trap around the MLMG solve (REEF3D's grid
// init has a benign NaN, so the global amrex.fpe_trap_invalid fires too early)
#include <AMReX_BLBackTrace.H>
#include <csignal>
#include <cfenv>

static void reef_set_fpe_trap(bool on)
{
#if defined(__APPLE__) && defined(__aarch64__)
    fenv_t env;
    fegetenv(&env);
    if(on) env.__fpcr |=  (__fpcr_trap_invalid | __fpcr_trap_divbyzero);
    else   env.__fpcr &= ~(__fpcr_trap_invalid | __fpcr_trap_divbyzero);
    fesetenv(&env);
#elif defined(__linux__)
    if(on) feenableexcept(FE_INVALID | FE_DIVBYZERO);
    else   fedisableexcept(FE_INVALID | FE_DIVBYZERO);
#else
    amrex::ignore_unused(on);
#endif
}

using namespace amrex;

amrex_solver::amrex_solver(lexer *p)
{
    // density selection mirrors pjm.cpp; the heat/concentration variants need
    // the heat*/concentration* pointers -> extend the constructor when wiring
    // those solver paths
    pd = nullptr;

    if(p->F80==0 && p->F300==0 && p->W90==0 && p->W30==0 && p->C10==0 && p->H10==0)
    {
        if(p->Q10==0)
        pd = new density_f(p);

        else if(p->Q10==1)
        pd = new density_pst(p);
    }

    if(p->H10==0 && p->W30==1 && p->F80==0 && p->F300==0 && p->W90==0)
    pd = new density_comp(p);

    if(p->F80>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
    pd = new density_vof(p);

    if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);

    if(p->F300>=1)
    pd = new density_rheo(p);

    if(pd==nullptr)
    pd = new density_f(p);

    gcval_press=40;

    gcval_u=7;
    gcval_v=8;
    gcval_w=9;
}

amrex_solver::~amrex_solver()
{
    delete pd;
}

// REEF_CF_FLUX: coarse/fine face-flux consistency probe.
//
// The composite projection is solvable only if the coarse face flux at every
// coarse-fine interface equals the area-weighted average of the fine faces
// covering it. fill_rhs builds rhs[lev] = -div(umac[lev]) from each level's own
// staged faces and nothing enforces that agreement, so any mismatch injects a
// source into the composite divergence that no pressure field can cancel: the
// solve then cannot converge and level 0's bottom solve fails every V-cycle.
//
// Only interface-NORMAL faces are compared -- a coarse face whose two adjacent
// coarse cells disagree on fine coverage. Faces on the physical domain boundary
// are skipped: a fine patch abutting a wall is a BC, not a C-F coupling. Faces
// shared by two coarse boxes are visited once per box, so `faces` is an upper
// bound on the count; the max-norms are unaffected.
static void reef_cf_flux_check(int nlev,
                               const Vector<BoxArray>& sgrids,
                               const Vector<Geometry>& sgeom,
                               const Vector<DistributionMapping>& dmaps,
                               const Vector<Array<MultiFab,AMREX_SPACEDIM>>& umac,
                               const IntVect& rr)
{
    // faces the fine patch does not cover keep this value, so an unfilled face
    // can never be mistaken for a perfectly matching one
    constexpr double sentinel = 1.0e30;

    for(int lev=1; lev<nlev; ++lev)
    {
        const BoxArray cba = amrex::coarsen(sgrids[lev], rr); // fine patch in coarse index space

        // ---- fine faces averaged down to coarse resolution ----
        Array<MultiFab,AMREX_SPACEDIM> favg;
        for(int d=0; d<AMREX_SPACEDIM; ++d)
        favg[d].define(amrex::convert(cba, IntVect::TheDimensionVector(d)), dmaps[lev], 1, 0);

        const Array<const MultiFab*,AMREX_SPACEDIM> fptr
            {AMREX_D_DECL(&umac[lev][0], &umac[lev][1], &umac[lev][2])};
        const Array<MultiFab*,AMREX_SPACEDIM> cptr
            {AMREX_D_DECL(&favg[0], &favg[1], &favg[2])};
        amrex::average_down_faces(fptr, cptr, rr, 0);

        // ---- restage them on the coarse level's own layout ----
        Array<MultiFab,AMREX_SPACEDIM> conf;
        for(int d=0; d<AMREX_SPACEDIM; ++d)
        {
            conf[d].define(amrex::convert(sgrids[lev-1], IntVect::TheDimensionVector(d)),
                           dmaps[lev-1], 1, 0);
            conf[d].setVal(sentinel);
            conf[d].ParallelCopy(favg[d], 0, 0, 1);
        }

        // ---- fine-coverage mask on the coarse level (1 ghost) ----
        iMultiFab cov(sgrids[lev-1], dmaps[lev-1], 1, 1);
        cov.setVal(0);
        {
            iMultiFab ones(cba, dmaps[lev], 1, 0);
            ones.setVal(1);
            cov.ParallelCopy(ones, 0, 0, 1);
        }
        cov.FillBoundary(sgeom[lev-1].periodicity());

        const Box cdom = sgeom[lev-1].Domain();

        for(int d=0; d<AMREX_SPACEDIM; ++d)
        {
            double dmax = 0.0, cmax = 0.0;
            long   nface = 0, nskip = 0;
            int    bi = 0, bj = 0, bk = 0;

            // conf/umac/cov share the coarse dmap and box order, so one MFIter
            // indexes all three (the face BoxArrays are convert()s of the cell one)
            for(MFIter mfi(conf[d]); mfi.isValid(); ++mfi)
            {
                const Box& fbx = mfi.validbox();
                auto const& cf = conf[d].const_array(mfi);
                auto const& uc = umac[lev-1][d].const_array(mfi);
                auto const& mk = cov.const_array(mfi);

                amrex::LoopOnCpu(fbx, [&] (int ii, int jj, int kk)
                {
                    const IntVect hi(AMREX_D_DECL(ii,jj,kk));
                    IntVect lo = hi; lo[d] -= 1;

                    if(hi[d] <= cdom.smallEnd(d) || hi[d] > cdom.bigEnd(d)) return; // domain face
                    if(mk(lo) == mk(hi))                                    return; // not a C-F face
                    if(cf(ii,jj,kk) >= 0.5*sentinel)         { ++nskip;      return; }

                    const double diff = std::abs(cf(ii,jj,kk) - uc(ii,jj,kk));
                    if(diff > dmax) { dmax = diff; bi = ii; bj = jj; bk = kk; }
                    cmax = std::max(cmax, std::abs(uc(ii,jj,kk)));
                    ++nface;
                });
            }

            const double local_dmax = dmax;
            double r[2] = {dmax, cmax};
            long   c[2] = {nface, nskip};
            ParallelDescriptor::ReduceRealMax(r,2);
            ParallelDescriptor::ReduceLongSum(c,2);

            amrex::Print()<<"CFFLUX lev "<<lev<<"->"<<lev-1<<" dir "<<d
                <<"  faces "<<c[0]<<"  uncovered "<<c[1]
                <<"  max|avg(fine)-coarse| "<<r[0]
                <<"  max|coarse| "<<r[1]<<std::endl;

            if(r[0] > 0.0 && local_dmax == r[0])
            {
                std::printf("CFFLUX   worst lev %d dir %d at (%d,%d,%d)\n", lev, d, bi, bj, bk);
                std::fflush(stdout);
            }
            ParallelDescriptor::Barrier();
        }
    }
}

void amrex_solver::setup(lexer *p, fdm *a, ghostcell *pgc, const field1 &u, const field2 &v, const field3 &w, const field4 &phi, double alpha)
{
    // Prerequisites: u/v/w and the level set must have current ghost cells
    // (pgc->start1/2/3/4 before calling) -- the staging below reads one ghost
    // layer for the low faces of each box, and roface reads the level set
    // through a->phi (the phi parameter documents the dependency).

    const int nlev = p->nlevs;

    // ---- solver-side hierarchy (see header: pseudo-2D y-doubling) ----
    // REEF_MLMG_NOYDOUBLE: A/B knob for the hidden-direction-at-nlev>1 path.
    // It does NOT currently work: the composite solve stalls (see the semi-
    // coarsening/hidden-direction comment in the LPInfo block below), so the
    // knob exists to re-test it after that is fixed, not as a production
    // setting. Forcing ydouble=false at nlev>1 aborts with
    // "MLMG: Failed to converge after 250 iterations, resid/resid0 = 3.4e-3".
    ydouble = (p->j_dir == 0) && nlev > 1;
    if(std::getenv("REEF_MLMG_NOYDOUBLE"))
    ydouble = false;

    sgeom.resize(nlev);
    sgrids.resize(nlev);
    sgeom[0]  = p->amrex_geometry[0];
    sgrids[0] = p->amrex_box_array[0];

    // Give the solver's BASE level two y-cells when y-doubling. AMReX's MG
    // coarsening test is Box::coarsenable(refrat, min_width), and its first
    // check is size().allGE(refrat*min_width) across ALL directions -- so a
    // 1-cell y makes the box un-coarsenable in x and z too, semicoarsening's
    // fallback then zeroes every ratio, and the hierarchy stops dead at one MG
    // level (measured: "# of MG levels on the coarsest AMR level: 1", with 95%
    // of the solve time in a bottom solve over the whole 192000-cell level 0).
    // setHiddenDirection is the only thing that zeroes min_width, and it is
    // unavailable here because it forces the AMR ratio anisotropic. Two y-cells
    // clear the min_width=2 bar; semicoarsening then pins y and lets x and z
    // coarsen normally.
    //
    // The extra plane is a copy, so the physical y-extent is DOUBLED rather than
    // the mesh refined: dy stays equal to dx and dz at every level, which keeps
    // the operator's cells isotropic and dhy comparable to dhx/dhz.
    if(ydouble)
    {
        const Geometry& g0 = p->amrex_geometry[0];

        Box d0 = g0.Domain();
        d0.setSmall(1, 0);
        d0.setBig(1, 1);

        RealBox rb0 = g0.ProbDomain();
        rb0.setHi(1, rb0.lo(1) + 2.0*g0.CellSize(1));

        const int isper[AMREX_SPACEDIM] =
            {AMREX_D_DECL(g0.isPeriodic(0), g0.isPeriodic(1), g0.isPeriodic(2))};
        sgeom[0].define(d0, &rb0, g0.Coord(), isper);

        BoxList bl0;
        for(int b=0; b<static_cast<int>(p->amrex_box_array[0].size()); ++b)
        {
            Box bb = p->amrex_box_array[0][b];
            bb.setSmall(1, 0);
            bb.setBig(1, 1);
            bl0.push_back(bb);
        }
        sgrids[0] = BoxArray(std::move(bl0));
    }

    // Solver-side refinement ratio -- geometry and box arrays MUST agree, because
    // MLLinOp takes no ref ratio: it infers one by coarsening the fine domain
    // until it matches the coarse one. A geometry refined 2x in a direction the
    // box arrays left at 1x never matches, the inference falls through and
    // yields ref ratio 8 -> mlmg_lin_cc_interp aborts.
    //
    // p->ref_vec is 1 in a thin direction, so a pseudo-2D hierarchy would be
    // refined ANISOTROPICALLY, (2,1,2), which is what setHiddenDirection needs.
    // That path is BROKEN: forcing ydouble=false at nlev>1 stalls the composite
    // solve at resid/resid0 ~ 3.4e-3 after 250 iterations (measured 2026-08-12,
    // first solve of step 1). rhs and resid0 match the working y-doubled run
    // exactly, so the staging and the RHS are fine and the operator is at fault.
    //
    // The cause is NOT the flux register, contrary to an earlier note here.
    // YAFluxRegister::FineAdd scales by dt/(fine_dx[d] * prod(rr)) and sums
    // prod(rr)/rr[d] fine faces; since fine_dx[d] = crse_dx[d]/rr[d], the total
    // is <f>*dt/crse_dx[d] for ANY rr, isotropic or not. The weighting is
    // dimensionally correct for (2,1,2) -- verified against
    // AMReX_YAFluxRegister.H ~542 and AMReX_YAFluxRegister_3D_K.H. MLCellLinOp
    // also passes AMRRefRatioVect (the IntVect, hidden dir = 1) to the register
    // at ~714, so it is not a scalar/vector mix-up either. The real cause is
    // still unidentified -- do not spend the search on reflux again.
    //
    // Refining the thin direction instead keeps the ratio isotropic and works.
    //
    // Refining the thin direction too keeps the ratio isotropic. The added
    // y-planes are exact copies of plane 0 (see the replication loop below) and
    // carry zero transverse velocity, so nothing physical changes; only the
    // solver-side index space grows.
    const IntVect rr = ydouble ? IntVect(p->ref_vec.max()) : p->ref_vec;

    for(int lev=1; lev<nlev; ++lev)
    {
        sgeom[lev] = amrex::refine(sgeom[lev-1], rr);

        if(ydouble)
        {
            // grid_amrex refined the box arrays in the active directions only and
            // left them one cell thick in y; stretch them to the y-extent that
            // sgeom[lev] now has, so grids and geometry agree. The base level
            // already carries two planes, hence lev+1.
            const int nyf = 1 << (lev+1);
            BoxList bl;
            for(int b=0; b<static_cast<int>(p->amrex_box_array[lev].size()); ++b)
            {
                Box bb = p->amrex_box_array[lev][b];
                bb.setSmall(1, 0);
                bb.setBig(1, nyf-1);
                bl.push_back(bb);
            }
            sgrids[lev] = BoxArray(std::move(bl));
        }
        else
        sgrids[lev] = p->amrex_box_array[lev];
    }

    // ---- (re)allocate the AMReX-side containers when the hierarchy changed ----
    bool rebuild = static_cast<int>(beta.size()) != nlev;
    for(int lev=0; lev<nlev && !rebuild; ++lev)
    rebuild = rhs[lev].boxArray() != sgrids[lev]
           || rhs[lev].DistributionMap() != p->amrex_distribution_mapping[lev];

    if(rebuild)
    {
        umac.clear();
        beta.clear();
        flux.clear();
        rhs.clear();
        pcorr.clear();
        umac.resize(nlev);
        beta.resize(nlev);
        flux.resize(nlev);
        rhs.resize(nlev);
        pcorr.resize(nlev);

        for(int lev=0; lev<nlev; ++lev)
        {
            for(int d=0; d<AMREX_SPACEDIM; ++d)
            {
                const BoxArray fba = convert(sgrids[lev], IntVect::TheDimensionVector(d));
                umac[lev][d].define(fba, p->amrex_distribution_mapping[lev], 1, 0);
                beta[lev][d].define(fba, p->amrex_distribution_mapping[lev], 1, 0);
                flux[lev][d].define(fba, p->amrex_distribution_mapping[lev], 1, 0);
            }
            rhs[lev].define(sgrids[lev], p->amrex_distribution_mapping[lev], 1, 0);
            pcorr[lev].define(sgrids[lev], p->amrex_distribution_mapping[lev], 1, 1);

            // warm-start seed: the solved potential is phi = alpha*dt*press
            // (see solve/pressure_update), so scale the stored physical
            // pressure by alpha*dt to seed it after an (re)allocation. Between
            // stages pcorr persists as phi and needs no reseeding.
            // press lives on the REEF3D (ny=1) grids; broadcast plane 0 into
            // all solver planes (same dmap and box order, so sharing the
            // MFIter is safe).
            const double adt = alpha*p->dt;
            pcorr[lev].setVal(0.0);
            for(MFIter mfi(pcorr[lev]); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                auto const& dst = pcorr[lev].array(mfi);
                auto const& src = a->press.GetMultiFab(lev).const_array(mfi);
                amrex::LoopOnCpu(vbx, [&] (int ii, int jj, int kk)
                {
                    dst(ii,jj,kk) = adt*src(ii, (ydouble ? 0 : jj), kk);
                });
            }
        }
    }

    // ---- stage face density and MAC velocity ----
    // REEF3D stores u/v/w cell-centred with the value living on the cell's
    // HIGH face; AMReX face index f is the LOW face of cell f, so face f in x
    // holds u(f-1). Every cell writes its three low faces; the last cell of a
    // tile row also writes its high face so the overlapping face planes of
    // adjacent boxes are filled in both FABs. The duplicate writes are
    // consistent because both evaluate the same exchanged cell data.
    // Faces touching a solid cell get beta = 0 (the no-flux fold) and a zero
    // velocity so their divergence contribution vanishes.
    MultiGridLOOP
    {
        auto const& bx = beta[p->level][0].atLocalIdx(p->amr_local_fab_idx).array();
        auto const& by = beta[p->level][1].atLocalIdx(p->amr_local_fab_idx).array();
        auto const& bz = beta[p->level][2].atLocalIdx(p->amr_local_fab_idx).array();
        auto const& uf = umac[p->level][0].atLocalIdx(p->amr_local_fab_idx).array();
        auto const& vf = umac[p->level][1].atLocalIdx(p->amr_local_fab_idx).array();
        auto const& wf = umac[p->level][2].atLocalIdx(p->amr_local_fab_idx).array();

        const int oi = p->amr_tile_lo.x;
        const int oj = p->amr_tile_lo.y;
        const int ok = p->amr_tile_lo.z;

        ILOOP
        JLOOP
        KLOOP
        {
            const int I = oi + i;
            const int J = oj + j;
            const int K = ok + k;

            const bool sc = p->flag4(i,j,k) <= SOLID_FLAG;

            // low faces
            if(sc || p->flag4(i-1,j,k) <= SOLID_FLAG)
            {
                bx(I,J,K) = 0.0;
                uf(I,J,K) = 0.0;
            }
            else
            {
                bx(I,J,K) = 1.0/pd->roface(p,a,-1,0,0);
                uf(I,J,K) = u(i-1,j,k);
            }

            if(ydouble)
            {
                // pseudo-2D: the solver planes are identical copies, so give
                // the y-faces the real cell density (the plane-to-plane
                // coupling lets the smoother damp plane-antisymmetric modes)
                // and zero transverse velocity
                by(I,J,K) = sc ? 0.0 : 1.0/pd->roface(p,a,0,0,0);
                vf(I,J,K) = 0.0;
            }
            else if(sc || p->flag4(i,j-1,k) <= SOLID_FLAG)
            {
                by(I,J,K) = 0.0;
                vf(I,J,K) = 0.0;
            }
            else
            {
                by(I,J,K) = 1.0/pd->roface(p,a,0,-1,0);
                vf(I,J,K) = v(i,j-1,k);
            }

            if(sc || p->flag4(i,j,k-1) <= SOLID_FLAG)
            {
                bz(I,J,K) = 0.0;
                wf(I,J,K) = 0.0;
            }
            else
            {
                bz(I,J,K) = 1.0/pd->roface(p,a,0,0,-1);
                wf(I,J,K) = w(i,j,k-1);
            }

            // high faces at the tile end close the face box
            if(i==IMAX_LOOP)
            {
                if(sc || p->flag4(i+1,j,k) <= SOLID_FLAG)
                {
                    bx(I+1,J,K) = 0.0;
                    uf(I+1,J,K) = 0.0;
                }
                else
                {
                    bx(I+1,J,K) = 1.0/pd->roface(p,a,1,0,0);
                    uf(I+1,J,K) = u(i,j,k);
                }
            }

            if(j==JMAX_LOOP)
            {
                if(ydouble)
                {
                    by(I,J+1,K) = by(I,J,K);
                    vf(I,J+1,K) = 0.0;
                }
                else if(sc || p->flag4(i,j+1,k) <= SOLID_FLAG)
                {
                    by(I,J+1,K) = 0.0;
                    vf(I,J+1,K) = 0.0;
                }
                else
                {
                    by(I,J+1,K) = 1.0/pd->roface(p,a,0,1,0);
                    vf(I,J+1,K) = v(i,j,k);
                }
            }

            if(k==KMAX_LOOP)
            {
                if(sc || p->flag4(i,j,k+1) <= SOLID_FLAG)
                {
                    bz(I,J,K+1) = 0.0;
                    wf(I,J,K+1) = 0.0;
                }
                else
                {
                    bz(I,J,K+1) = 1.0/pd->roface(p,a,0,0,1);
                    wf(I,J,K+1) = w(i,j,k);
                }
            }
        }
    }

    // ---- replicate plane 0 across the added y-planes ----
    // level 0 now carries two planes as well, so this starts at 0
    if(ydouble)
    for(int lev=0; lev<nlev; ++lev)
    for(int d=0; d<AMREX_SPACEDIM; ++d)
    for(MFIter mfi(beta[lev][d]); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& b = beta[lev][d].array(mfi);
        auto const& u = umac[lev][d].array(mfi);
        amrex::LoopOnCpu(vbx, [&] (int ii, int jj, int kk)
        {
            if(jj > 0)
            {
                b(ii,jj,kk) = b(ii,0,kk);
                u(ii,jj,kk) = u(ii,0,kk);
            }
        });
    }

    // REEF_CF_FLUX=1: does the coarse face flux match the average of the fine
    // faces at every C-F interface? A mismatch makes the composite rhs
    // unsolvable -- see reef_cf_flux_check.
    if(nlev > 1 && std::getenv("REEF_CF_FLUX"))
    reef_cf_flux_check(nlev, sgrids, sgeom, p->amrex_distribution_mapping, umac, rr);

    // ---- build the composite operator and the multigrid solver ----
    // The operator structure (grids, geometry, BC types, hidden direction) and
    // the MLMG object depend only on the grid hierarchy, so they are rebuilt
    // only on the first call and on regrid -- not on every RK stage. This
    // mirrors hypre_ssamg's solver_created/created_nlevs guard. `rebuild` is
    // already true whenever the MultiFab containers above were reallocated
    // (grid or dmap change); the extra flags make the intent explicit and
    // catch a stale linop even if the containers happened to match.
    // REEF_MLMG_FORCE_REBUILD: A/B knob -- rebuild the operator every stage as
    // before the guard, to measure the setup overhead the guard removes.
    const bool solver_rebuild = rebuild || !solver_created || created_nlevs != nlev
                              || std::getenv("REEF_MLMG_FORCE_REBUILD");

    if(solver_rebuild)
    {
        Vector<DistributionMapping> dmaps(nlev);
        for(int lev=0; lev<nlev; ++lev)
        dmaps[lev] = p->amrex_distribution_mapping[lev];

        LPInfo info;
        info.setAgglomeration(true);
        info.setConsolidation(true);
        // REEF3D runs pseudo-2D cases with a single cell in y. AMReX coarsens
        // isotropically, so a 1-cell direction blocks coarsening in every
        // direction: the hierarchy collapses to a single MG level and each MLMG
        // iteration degenerates into an unpreconditioned BiCGStab over the whole
        // domain. There are two ways out and they are mutually exclusive
        // (MLLinOp asserts), so exactly one is enabled here:
        //
        //  - hidden direction: sets the thin axis' coarsen ratio to 1 and its min
        //    width to 0. Used at nlev==1, where it converges in a handful of
        //    V-cycles. NOT usable with AMR levels: it forces the AMR ratio
        //    anisotropic, which reflux mis-weights (see the rr comment above),
        //    and it is what the three nlev>1-only local AMReX patches exist for.
        //  - semicoarsening: coarsens the active directions once the thin one
        //    bottoms out. Needs BOTH calls -- MLLinOp gates every step on
        //    numsclevs < max_semicoarsening_level, which defaults to 0, so
        //    setSemicoarsening alone is silently a no-op.
        //
        // With AMR levels the grids are y-doubled to keep the ratio isotropic, so
        // the hidden direction is off and semicoarsening carries the hierarchy.
        const int ny = sgeom[0].Domain().length(1);
        if(ny==1 && !ydouble && !std::getenv("REEF_MLMG_NOHIDDEN"))
        {
            info.setHiddenDirection(1);
        }
        else
        {
            const char* msc = std::getenv("REEF_MLMG_SEMICOARSEN");
            info.setSemicoarsening(true);
            info.setMaxSemicoarseningLevel(msc ? std::atoi(msc) : 8);
        }

        // bisection knobs for multi-level debugging
        if(std::getenv("REEF_MLMG_NOAGG"))
        {
            info.setAgglomeration(false);
            info.setConsolidation(false);
        }
        if(const char* mc = std::getenv("REEF_MLMG_MAXCOARSE"))
        info.setMaxCoarseningLevel(std::atoi(mc));

        linop = std::make_unique<MLABecLaplacian>(sgeom, sgrids, dmaps, info);

        // REEF3D pressure BC (gcval 40) is Neumann on every physical boundary;
        // the all-Neumann nullspace is handled inside MLMG
        linop->setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann, LinOpBCType::Neumann, LinOpBCType::Neumann)},
                           {AMREX_D_DECL(LinOpBCType::Neumann, LinOpBCType::Neumann, LinOpBCType::Neumann)});

        for(int lev=0; lev<nlev; ++lev)
        linop->setLevelBC(lev, nullptr);

        // operator is -div(beta grad) with beta = 1/ro on faces; dt stays in the
        // RHS and in the velocity correction, matching the pjm/pjm_corr convention
        linop->setScalars(0.0, 1.0);

        for(int lev=0; lev<nlev; ++lev)
        linop->setACoeffs(lev, 0.0);

        mlmg = std::make_unique<MLMG>(*linop);
        mlmg->setMaxIter(p->N46);
        // tolerance p->N44 is passed at solve time

        if(const char* mv = std::getenv("REEF_MLMG_VERBOSE"))
        {
            mlmg->setVerbose(std::atoi(mv));
            mlmg->setBottomVerbose(std::atoi(mv));
        }
        if(const char* fi = std::getenv("REEF_MLMG_FIXEDITER"))
        mlmg->setFixedIter(std::atoi(fi));
        if(std::getenv("REEF_MLMG_BOTTOM_SMOOTHER"))
        mlmg->setBottomSolver(MLMG::BottomSolver::smoother);

        // Bottom solver: BoomerAMG via hypre. AMReX' default bicgstab breaks down
        // on the singular all-Neumann bottom problem, and the semicoarsened
        // hierarchy leaves a bottom grid whose x:z aspect ratio mirrors the
        // domain's -- exactly the anisotropy an algebraic coarsening handles and
        // a geometric one does not. Set before the env block so REEF_MLMG_BOTTOM
        // can still override it for bisection.
        mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
        mlmg->setBottomTolerance(1.e-3);   // inside a V-cycle; AMReX default 1e-4
        mlmg->setBottomMaxIter(100);

        // REEF_MLMG_BOTTOM=cg|bicgstab|hypre|smoother
        if(const char* bs = std::getenv("REEF_MLMG_BOTTOM"))
        {
            const std::string s(bs);
            if(s=="cg")            mlmg->setBottomSolver(MLMG::BottomSolver::cg);
            else if(s=="bicgstab") mlmg->setBottomSolver(MLMG::BottomSolver::bicgstab);
            else if(s=="smoother") mlmg->setBottomSolver(MLMG::BottomSolver::smoother);
            else if(s=="hypre")    mlmg->setBottomSolver(MLMG::BottomSolver::hypre);
        }
        if(const char* bi = std::getenv("REEF_MLMG_BOTTOM_MAXITER"))
        mlmg->setBottomMaxIter(std::atoi(bi));

        if(std::getenv("REEF_MLMG_BOTTOM_SIZE"))
        {
            mlmg->setNSolve(1);              // bottom problem on a single grid/rank
            mlmg->setNSolveGridSize(16);
            mlmg->setBottomTolerance(1.e-3);
        }

        // F-cycles for the first iterations. AMReX defaults max_fmg_iters to 0,
        // i.e. plain V-cycles with no F-cycle start-up; REEF_MLMG_MAXFMGITER=0
        // restores that for comparison.
        {
            const char* fm = std::getenv("REEF_MLMG_MAXFMGITER");
            mlmg->setMaxFmgIter(fm ? std::atoi(fm) : 4);
        }

        // Heavier smoothing than AMReX' default nu1=nu2=2. Each V-cycle costs
        // more, but the density jump at the free surface makes the iteration
        // count the dominant term. REEF_MLMG_SMOOTH=2 restores the default.
        {
            const char* sm = std::getenv("REEF_MLMG_SMOOTH");
            const int nu = sm ? std::atoi(sm) : 3;
            mlmg->setPreSmooth(nu);
            mlmg->setPostSmooth(nu);
        }

        solver_created = true;
        created_nlevs  = nlev;
    }

    // Face coefficients change every solve with the moving free surface, so
    // rebind them on the (possibly reused) operator. MLMG::solve re-averages
    // the updated coefficients down the hierarchy via linop->prepareForSolve.
    for(int lev=0; lev<nlev; ++lev)
    linop->setBCoeffs(lev, GetArrOfConstPtrs(beta[lev]));

    // REEF_MLMG_APPLY: operator sanity probe -- L(1) must be ~0 for the
    // all-Neumann operator (zero row sums); NaN here means the operator/BC
    // setup is broken, clean means the smoother path is at fault
    if(std::getenv("REEF_MLMG_APPLY"))
    {
        Vector<MultiFab> xin(nlev), xout(nlev);
        for(int lev=0; lev<nlev; ++lev)
        {
            xin[lev].define(sgrids[lev], p->amrex_distribution_mapping[lev], 1, 1);
            xout[lev].define(sgrids[lev], p->amrex_distribution_mapping[lev], 1, 0);
            xin[lev].setVal(1.0);
            xout[lev].setVal(0.0);
        }
        MLMG probe(*linop);
        probe.apply(GetVecOfPtrs(xout), GetVecOfPtrs(xin));
        for(int lev=0; lev<nlev; ++lev)
        amrex::Print()<<"CHECK apply L(1) lev "<<lev
            <<"  ["<<xout[lev].min(0)<<","<<xout[lev].max(0)<<"] nan="<<xout[lev].contains_nan()
            <<std::endl;
    }

    // REEF_MLMG_CHECK: staged-input diagnostics (min/max/NaN per level)
    if(std::getenv("REEF_MLMG_CHECK"))
    {
        // adaptive-regrid boxes (blocking_factor 1) may not be 2-coarsenable;
        // in release builds MLLinOp's check is a no-op assert, so a violation
        // corrupts the C-F machinery silently
        for(int lev=1; lev<nlev; ++lev)
        {
            amrex::Print()<<"CHECK lev "<<lev<<" boxes "<<sgrids[lev].size()
                <<" coarsenable(2)="<<sgrids[lev].coarsenable(2)
                <<" min_width="<<sgrids[lev].minimalBox().shortside()
                <<std::endl;
            for(int b=0; b<static_cast<int>(sgrids[lev].size()); ++b)
            amrex::Print()<<"CHECK lev "<<lev<<" box "<<b<<" "<<sgrids[lev][b]<<std::endl;
        }
        // zero-diagonal probe: cells whose every face has beta=0 divide the
        // GSRB smoother by zero (0/0 -> NaN, spreads through the zero faces)
        for(int lev=0; lev<nlev; ++lev)
        {
            long nzero = 0;
            for(MFIter mfi(rhs[lev]); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                auto const& bx_ = beta[lev][0].const_array(mfi);
                auto const& by_ = beta[lev][1].const_array(mfi);
                auto const& bz_ = beta[lev][2].const_array(mfi);
                amrex::LoopOnCpu(vbx, [&] (int ii, int jj, int kk)
                {
                    const double s = bx_(ii,jj,kk)+bx_(ii+1,jj,kk)
                                   + by_(ii,jj,kk)+by_(ii,jj+1,kk)
                                   + bz_(ii,jj,kk)+bz_(ii,jj,kk+1);
                    if(s==0.0) ++nzero;
                });
            }
            amrex::ParallelDescriptor::ReduceLongSum(nzero);
            amrex::Print()<<"CHECK zero-diagonal cells lev "<<lev<<": "<<nzero<<std::endl;
        }
        // flag4 census: the no-flux fold keys off flag4<=SOLID_FLAG, so a solid
        // region carrying any other flag is silently solved as fluid
        {
            long nsol=0, nzero_flag=0, nair=0, nwat=0, nio=0, nother=0;
            MultiGridLOOP
            {
                ILOOP JLOOP KLOOP
                {
                    const int f = p->flag4(i,j,k);
                    if(f<=SOLID_FLAG)                        ++nsol;
                    else if(f==0)                            ++nzero_flag;
                    else if(f==AIR_FLAG)                     ++nair;
                    else if(f>=WATER_FLAG)                   ++nwat;
                    else if(f==INFLOW_FLAG||f==OUTFLOW_FLAG) ++nio;
                    else                                     ++nother;
                }
            }
            long c[6] = {nsol,nzero_flag,nair,nwat,nio,nother};
            amrex::ParallelDescriptor::ReduceLongSum(c,6);
            amrex::Print()<<"CHECK flag4 census  solid(<=-19)="<<c[0]<<"  zero="<<c[1]
                <<"  air="<<c[2]<<"  water="<<c[3]<<"  in/out="<<c[4]<<"  other="<<c[5]<<std::endl;
        }
        for(int lev=0; lev<nlev; ++lev)
        for(int d=0; d<AMREX_SPACEDIM; ++d)
        {
            amrex::Print()<<"CHECK setup lev "<<lev<<" dir "<<d
                <<"  beta["<<beta[lev][d].min(0)<<","<<beta[lev][d].max(0)<<"] nan="<<beta[lev][d].contains_nan()
                <<"  umac["<<umac[lev][d].min(0)<<","<<umac[lev][d].max(0)<<"] nan="<<umac[lev][d].contains_nan()
                <<std::endl;
        }
        for(int lev=0; lev<nlev; ++lev)
        amrex::Print()<<"CHECK setup lev "<<lev
            <<"  pcorr["<<pcorr[lev].min(0)<<","<<pcorr[lev].max(0)<<"] nan="<<pcorr[lev].contains_nan(0,1,1)
            <<std::endl;
    }
}

void amrex_solver::fill_rhs(lexer *p)
{
    // The solve variable is the projection potential phi = alpha*dt*p, which
    // satisfies -div(beta grad phi) = -div(umac). Keeping alpha*dt out of the
    // RHS (and out of the solution) makes phi independent of the RK stage, so
    // the carried-over pcorr is a genuine warm start and bnorm no longer swings
    // with alpha. The alpha*dt factor is restored in pressure_update. Solid
    // cells get rhs = 0 automatically because all their faces were zeroed.
    for(int lev=0; lev<p->nlevs; ++lev)
    {
        const auto dxinv = p->amrex_geometry[lev].InvCellSizeArray();
        const Real dxi = dxinv[0];
        const Real dyi = dxinv[1];
        const Real dzi = dxinv[2];

        for(MFIter mfi(rhs[lev],MFIter_TILING); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            auto const& r  = rhs[lev].array(mfi);
            auto const& uf = umac[lev][0].const_array(mfi);
            auto const& vf = umac[lev][1].const_array(mfi);
            auto const& wf = umac[lev][2].const_array(mfi);

            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int ii, int jj, int kk) noexcept
            {
                r(ii,jj,kk) = -( (uf(ii+1,jj,kk)-uf(ii,jj,kk))*dxi
                               + (vf(ii,jj+1,kk)-vf(ii,jj,kk))*dyi
                               + (wf(ii,jj,kk+1)-wf(ii,jj,kk))*dzi );
            });
        }
    }
}

double amrex_solver::solve(lexer *p)
{
    if(std::getenv("REEF_MLMG_CHECK"))
    for(int lev=0; lev<p->nlevs; ++lev)
    amrex::Print()<<"CHECK rhs lev "<<lev
        <<"  ["<<rhs[lev].min(0)<<","<<rhs[lev].max(0)<<"] nan="<<rhs[lev].contains_nan()
        <<std::endl;

    const bool fpe_trap = std::getenv("REEF_MLMG_FPE") != nullptr;
    if(fpe_trap)
    {
        std::signal(SIGFPE, amrex::BLBackTrace::handler);
        std::signal(SIGILL, amrex::BLBackTrace::handler); // Apple Silicon raises SIGILL for FP traps
        reef_set_fpe_trap(true);
    }

    const double finres = mlmg->solve(GetVecOfPtrs(pcorr), GetVecOfConstPtrs(rhs), p->N44, 0.0);

    if(fpe_trap)
    reef_set_fpe_trap(false);

    p->solveriter = mlmg->getNumIters();

    if(std::getenv("REEF_MLMG_CHECK"))
    for(int lev=0; lev<p->nlevs; ++lev)
    amrex::Print()<<"CHECK sol lev "<<lev
        <<"  pcorr["<<pcorr[lev].min(0)<<","<<pcorr[lev].max(0)<<"] nan="<<pcorr[lev].contains_nan(0,1,1)
        <<std::endl;

    // REEF_MLMG_WALLROW: print the level-0 solution along the bottom row
    if(std::getenv("REEF_MLMG_WALLROW"))
    {
        for(MFIter mfi(pcorr[0]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& s = pcorr[0].const_array(mfi);
            for(int ii=20; ii<=45; ++ii)
            if(vbx.contains(IntVect(ii,0,0)))
            std::printf("WALLROW p(%d,0,0)=%.6e  p(%d,0,1)=%.6e\n",
                        ii, s(ii,0,0), ii, s(ii,0,1));
            std::fflush(stdout);
        }
        if(p->nlevs > 1)
        for(MFIter mfi(beta[1][0]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& b = beta[1][0].const_array(mfi);
            for(int ii=48; ii<=71; ii+=2)
            if(vbx.contains(IntVect(ii,0,0)))
            std::printf("WALLROW betax1(%d,0,0)=%.6e  betax1(%d,0,40)=%.6e\n",
                        ii, b(ii,0,0), ii, vbx.contains(IntVect(ii,0,40)) ? b(ii,0,40) : -1.0);
            std::fflush(stdout);
        }
    }

    // REEF_MLMG_RESLOC: locate the worst composite-residual cells after the solve
    if(std::getenv("REEF_MLMG_RESLOC"))
    {
        Vector<MultiFab> Lp(p->nlevs);
        for(int lev=0; lev<p->nlevs; ++lev)
        {
            Lp[lev].define(sgrids[lev], p->amrex_distribution_mapping[lev], 1, 0);
            Lp[lev].setVal(0.0);
        }
        MLMG probe(*linop);
        probe.apply(GetVecOfPtrs(Lp), GetVecOfPtrs(pcorr));

        for(int lev=0; lev<p->nlevs; ++lev)
        {
            double rmax = 0.0;
            int ri=-1, rj=-1, rk=-1;
            for(MFIter mfi(rhs[lev]); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();
                auto const& rr = rhs[lev].const_array(mfi);
                auto const& ll = Lp[lev].const_array(mfi);
                amrex::LoopOnCpu(vbx, [&] (int ii, int jj, int kk)
                {
                    const double r = std::fabs(rr(ii,jj,kk) - ll(ii,jj,kk));
                    if(r > rmax) { rmax = r; ri=ii; rj=jj; rk=kk; }
                });
            }
            std::printf("RESLOC rank %d lev %d max|r| %.4e at (%d,%d,%d)\n",
                        amrex::ParallelDescriptor::MyProc(), lev, rmax, ri, rj, rk);
            std::fflush(stdout);
        }
    }

    return finres;
}

void amrex_solver::ucorr(lexer *p, fdm *a, ghostcell *pgc, field1 &u, field2 &v, field3 &w)
{
    // flux = -(1/ro)grad(pcorr) evaluated with the operator's own discrete
    // gradient (including the asymmetric C-F stencils and the beta=0 solid
    // faces) -- the exact adjoint of the divergence the solve saw, so the
    // corrected field is discretely divergence-free on the composite grid.
    // pcorr is the potential phi = alpha*dt*p, so this flux already carries the
    // alpha*dt factor: the correction is u += flux, with no further scaling.
    // Zero first: with a hidden (pseudo-2D) direction MLMG does not write
    // that direction's flux, and the correction below must add 0 there.
    for(int lev=0; lev<p->nlevs; ++lev)
    for(int d=0; d<AMREX_SPACEDIM; ++d)
    flux[lev][d].setVal(0.0);

    mlmg->getFluxes(GetVecOfArrOfPtrs(flux), MLMG::Location::FaceCenter);

    if(std::getenv("REEF_MLMG_CHECK"))
    for(int lev=0; lev<p->nlevs; ++lev)
    for(int d=0; d<AMREX_SPACEDIM; ++d)
    amrex::Print()<<"CHECK flux lev "<<lev<<" dir "<<d
        <<"  ["<<flux[lev][d].min(0)<<","<<flux[lev][d].max(0)<<"] nan="<<flux[lev][d].contains_nan()
        <<std::endl;

    // REEF_MLMG_CHECK: track the largest applied velocity correction
    const bool chk = std::getenv("REEF_MLMG_CHECK") != nullptr;
    double cmax = 0.0;
    int cl=-1, cd=-1, ci=-1, cj=-1, ck=-1;

    // cell (i,j,k) stores its HIGH face -> AMReX face index (I+1,J,K) etc.
    // IULOOP/JVLOOP/KWLOOP exclude the domain high boundary face (BC-owned),
    // and the flag checks skip solid/ghost faces, mirroring pjm's u/v/wcorr
    MultiGridLOOP
    {
        auto const& fx = flux[p->level][0].atLocalIdx(p->amr_local_fab_idx).const_array();
        auto const& fy = flux[p->level][1].atLocalIdx(p->amr_local_fab_idx).const_array();
        auto const& fz = flux[p->level][2].atLocalIdx(p->amr_local_fab_idx).const_array();

        const int oi = p->amr_tile_lo.x;
        const int oj = p->amr_tile_lo.y;
        const int ok = p->amr_tile_lo.z;

        IULOOP
        JLOOP
        KLOOP
        UCHECK
        {
            const double du = fx(oi+i+1, oj+j, ok+k);
            u(i,j,k) += du;
            if(chk && std::fabs(du) > cmax)
            { cmax=std::fabs(du); cl=p->level; cd=0; ci=oi+i; cj=oj+j; ck=ok+k; }
        }

        ILOOP
        JVLOOP
        KLOOP
        VCHECK
        v(i,j,k) += fy(oi+i, oj+j+1, ok+k);

        ILOOP
        JLOOP
        KWLOOP
        WCHECK
        {
            const double dw = fz(oi+i, oj+j, ok+k+1);
            w(i,j,k) += dw;
            if(chk && std::fabs(dw) > cmax)
            { cmax=std::fabs(dw); cl=p->level; cd=2; ci=oi+i; cj=oj+j; ck=ok+k; }
        }
    }

    if(chk)
    {
        std::printf("UCORR rank %d max|du| %.4e lev %d dir %d cell (%d,%d,%d)\n",
                    amrex::ParallelDescriptor::MyProc(), cmax, cl, cd, ci, cj, ck);
        std::fflush(stdout);
    }

    pgc->start1(p,u,gcval_u);
    pgc->start2(p,v,gcval_v);
    pgc->start3(p,w,gcval_w);

    // the projection guarantees C-F interface flux consistency; covered
    // coarse faces still need the fine average
    for(int lev=p->nlevs-2; lev>=0; --lev)
    {
        u.average_down_level(p,lev);
        v.average_down_level(p,lev);
        w.average_down_level(p,lev);
    }
}

void amrex_solver::pressure_update(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
    // non-incremental solve: pcorr is the potential phi = alpha*dt*p, so the
    // physical pressure is phi/(alpha*dt) (all-Neumann solves are pinned to
    // zero mean inside MLMG; add a reference-pressure shift here if a
    // pressure_reference-style gauge is needed). press lives on the REEF3D
    // grids; with y-doubling active its box is the j=0 plane of the solver
    // box, so an index-wise copy reads back exactly that plane.
    const double iadt = 1.0/(alpha*p->dt);
    for(int lev=0; lev<p->nlevs; ++lev)
    {
        auto& pr = a->press.GetMultiFab(lev);
        for(MFIter mfi(pr); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const& dst = pr.array(mfi);
            auto const& src = pcorr[lev].const_array(mfi);
            amrex::LoopOnCpu(vbx, [&] (int ii, int jj, int kk)
            {
                dst(ii,jj,kk) = iadt*src(ii,jj,kk);
            });
        }
    }

    pgc->start4(p,a->press,gcval_press);
}

void amrex_solver::start(lexer *p, fdm *a, ghostcell *pgc, field1 &u, field2 &v, field3 &w, const field4 &phi, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    std::cout<<".";

    const double total_start = pgc->timer();

    // mirror pjm::vel_setup: the caller updates the velocity ghosts only
    // after the projection, but the staging reads one ghost layer
    double block_start = pgc->timer();
    pgc->start1(p,u,gcval_u);
    pgc->start2(p,v,gcval_v);
    pgc->start3(p,w,gcval_w);
    const double gc_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    setup(p,a,pgc,u,v,w,phi,alpha);
    const double setup_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    fill_rhs(p);
    const double fill_rhs_time = pgc->timer() - block_start;

    const double starttime = pgc->timer();

    solve(p);

    const double endtime = pgc->timer();

    block_start = pgc->timer();
    pressure_update(p,a,pgc,alpha);
    const double pressure_update_time = pgc->timer() - block_start;

    block_start = pgc->timer();
    ucorr(p,a,pgc,u,v,w);
    const double ucorr_time = pgc->timer() - block_start;

    p->poissoniter = p->solveriter;
    p->poissontime = endtime-starttime;

    const double total_time = pgc->timer() - total_start;

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
        std::cout<<"piter: "<<p->solveriter<<"  ptime: "<<std::setprecision(3)<<p->poissontime<<std::endl;

        if(std::getenv("REEF_MLMG_TIMING"))
        {
            const double denom = (total_time > 0.0) ? total_time : 1.0;
            std::cout<<"amrex_solver runtime breakdown (s | %total):"<<std::endl;
            std::cout<<std::setprecision(6)
                     <<"  gc: "<<gc_time<<" | "<<(100.0*gc_time/denom)<<"%\n"
                     <<"  setup: "<<setup_time<<" | "<<(100.0*setup_time/denom)<<"%\n"
                     <<"  fill_rhs: "<<fill_rhs_time<<" | "<<(100.0*fill_rhs_time/denom)<<"%\n"
                     <<"  solve: "<<p->poissontime<<" | "<<(100.0*p->poissontime/denom)<<"%\n"
                     <<"  pressure_update: "<<pressure_update_time<<" | "<<(100.0*pressure_update_time/denom)<<"%\n"
                     <<"  ucorr: "<<ucorr_time<<" | "<<(100.0*ucorr_time/denom)<<"%\n"
                     <<"  total: "<<total_time<<std::endl;
        }
    }
}
#endif
