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

#include"density.h"
#include"lexer.h"
#include"fdm.h"
#include"looping.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#if USE_AMREX
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#endif

#if USE_AMREX
// Make the face density single-valued across every coarse-fine interface, mirroring what
// AMReX does for the MLABecLaplacian b-coefficients: MLABecLaplacianT::averageDownCoeffs()
// -> averageDownCoeffsToCoarseAmrLevel() -> amrex::average_down_faces(). Because that call
// builds its temporary on coarsen(fine.boxArray()) in FACE index type, the coarse faces it
// overwrites include the C-F interface face itself, not just the faces interior to the
// covered region. That is the property REEF3D was missing: with roface() recomputed from phi
// independently on each level, a smoothed-Heaviside band crossing the C-F gave the coarse and
// fine sides of one physical face two different densities.
//
// AMReX averages the DIFFUSION COEFFICIENT beta = 1/rho, and mean(1/rho) != 1/mean(rho). We
// store rho (so every consumer stays a plain division and single-level results are unchanged),
// so the reciprocal is taken on the way into the temporaries and undone on the way back: the
// stored coarse value becomes 1/mean(1/rho_fine), and the consumer's 1/rho_stored is then
// exactly AMReX's mean(beta).
//
// The coarse temporary is seeded with a 0.0 sentinel rather than the coarse rho, so the
// write-back can touch only the faces average_down_faces actually filled. Seeding with rho
// and writing everything back would round-trip every uncovered coarse face through two
// reciprocals and perturb it by an ulp.
static void cf_sync_faces(lexer* p, fdm* a)
{
    if(p->nlevs <= 1) return;

    field* rof[AMREX_SPACEDIM] = {&a->rofx, &a->rofy, &a->rofz};

    for(int lev = p->nlevs-1; lev >= 1; --lev)
    {
        const int clev = lev-1;
        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> ffine, fcrse;

        for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

            ffine[dir].define(amrex::convert(p->amrex_box_array[lev],  e),
                              p->amrex_distribution_mapping[lev],  1, 0);
            fcrse[dir].define(amrex::convert(p->amrex_box_array[clev], e),
                              p->amrex_distribution_mapping[clev], 1, 0);

            // cell-stored HIGH face -> face MultiFab: face f lies between cells f-e and f,
            // i.e. it is the high face of cell f-e. Inverted to beta = 1/rho for the average.
            auto& r_mf = rof[dir]->GetMultiFab(lev);
            for(amrex::MFIter mfi(ffine[dir]); mfi.isValid(); ++mfi)
            {
                const amrex::Box& fbx = mfi.validbox();
                auto       f  = ffine[dir].array(mfi);
                const auto rr = r_mf.const_array(mfi);
                amrex::LoopOnCpu(fbx, [&] (int ii, int jj, int kk) noexcept
                { f(ii,jj,kk) = 1.0/rr(ii-e[0], jj-e[1], kk-e[2]); });
            }

            fcrse[dir].setVal(0.0);   // sentinel: untouched faces stay 0
        }

        const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> ffp{&ffine[0],&ffine[1],&ffine[2]};
        const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>       fcp{&fcrse[0],&fcrse[1],&fcrse[2]};

        amrex::average_down_faces(ffp, fcp, p->ref_vec, p->amrex_geometry[clev]);

        // REEF_CF_DENSITY_PROBE: max |rho_coarse - 1/mean(1/rho_fine)| over the covered faces,
        // i.e. exactly the two-valued-density gap this sync closes. Reported BEFORE the
        // write-back, since afterwards it is 0 by construction. Note levelset_RK3::phi_cf_probe
        // measures the same physical quantity but recomputes it from phi, so it keeps reporting
        // the raw gap and is unaffected by this fix -- use this number, not that one, to confirm
        // the sync fired. W1-W3 scale here means the Heaviside band is straddling the interface.
        const bool dprobe = (std::getenv("REEF_CF_DENSITY_PROBE") != nullptr);
        double gap[AMREX_SPACEDIM] = {0.0, 0.0, 0.0};

        for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if(dir==1 && p->j_dir!=1) continue;

            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
            auto& r_mf = rof[dir]->GetMultiFab(clev);
            for(amrex::MFIter mfi(r_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                auto       rr = r_mf.array(mfi);
                const auto f  = fcrse[dir].const_array(mfi);
                amrex::LoopOnCpu(bx, [&] (int ii, int jj, int kk) noexcept
                {
                    const double b = f(ii+e[0], jj+e[1], kk+e[2]);
                    if(b <= 0.0) return;                // uncovered face: keep the coarse value
                    const double rfine = 1.0/b;
                    if(dprobe) gap[dir] = std::max(gap[dir], std::fabs(rr(ii,jj,kk) - rfine));
                    rr(ii,jj,kk) = rfine;               // covered face: take the fine average
                });
            }
        }

        if(dprobe)
        {
            amrex::ParallelDescriptor::ReduceRealMax(gap, AMREX_SPACEDIM);
            amrex::Print()<<"  [cfdensity] clev="<<clev
                          <<"  max|rho_coarse - fine_avg| over covered faces  x/y/z = "
                          <<gap[0]<<" / "<<gap[1]<<" / "<<gap[2]
                          <<"   (W1-W3 = "<<(p->W1-p->W3)<<")"<<std::endl;
        }
    }
}
#endif

void density::update_faces(lexer* p, fdm* a)
{
    // Range -1..IMAX, NOT ULOOP/VLOOP/WLOOP. The staggered loops carry flag checks that skip
    // solid-adjacent faces and the domain-high face, but poisson_pcorr evaluates the density on
    // all six faces of every fluid cell -- including the ones it then folds into the diagonal.
    // An unfilled 0.0 there divides to inf and the fold turns it into NaN. The low ghost layer
    // is written for a different reason: a cell's LOW face is stored on its i-1 neighbour, which
    // at a tile/domain edge is a ghost cell. Nothing reads above IMAX, so the high side stops
    // there -- BLOOP's extra layer would only read phi two ghosts deep for no consumer.
    //
    // MFIter tiling is live (definitions_amrex.h: TilingIfNotGPU), so an interior tile's -1 layer
    // lands on the previous tile's valid cell in the same FAB. That is a redundant write, not a
    // corruption: both tiles evaluate the same phi pair for that face and get the same number.
    //
    // No start1/2/3 afterwards. Every value derives from phi, whose ghosts (domain BC, box halo,
    // and the coarse-fine ring) the caller has already filled, so the halo is right by
    // construction. A Neumann fill would in fact be WRONG at the domain-boundary low faces: it
    // would replace roface(phi(-1),phi(0)) with roface(phi(0),phi(1)), and that value multiplies
    // the Dirichlet pressure move in poisson_pcorr.
    MultiGridLOOP
    for(i = -1; i <= IMAX_LOOP; ++i)
    for(j = -1; j <= JMAX_LOOP; ++j)
    for(k = -1; k <= KMAX_LOOP; ++k)
    {
        a->rofx(i,j,k) = roface(p,a, 1,0,0);
        a->rofy(i,j,k) = roface(p,a, 0,1,0);
        a->rofz(i,j,k) = roface(p,a, 0,0,1);
    }

#if USE_AMREX
    cf_sync_faces(p,a);
#endif
}
