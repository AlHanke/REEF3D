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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"heaviside_ls.h"
#if USE_AMREX
#include "density_f.h"
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <cstdlib>
#include <iostream>
#endif

void initialize::hydrostatic(lexer *p, fdm *a, ghostcell *pgc)
{
    double maxh=0.0;
    BASELOOP
    maxh=MAX(maxh, p->pos_z());

    maxh=pgc->globalmax(maxh);
	
    if(p->F30>0)
    maxh=p->phimean;
	
	if(p->I12==1 && (p->I30==0||p->B90==0))
    {
        BASELOOP
        a->press(i,j,k) = (p->phimean-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
    }
	else if(p->I12==2 && (p->I30==0||p->B90==0))
    {
        BASELOOP
        a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22);
    }
	else if(p->I12==3 && (p->I30==0||p->B90==0))
    {
        BASELOOP
        a->press(i,j,k) = (maxh-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
    }
    else if(p->I12==4 && p->Y9==1)
    {
        #if USE_AMREX
        const double psi = p->psi;
        auto dens = density_f(p);

        // Column-by-column hydrostatic build. Each (i,j) column is anchored at its tile-bottom
        // cell with the analytic absolute value, then integrated UPWARD using the SAME face
        // density momentum/wpgrad use (dens.roface, field phi). That makes grad(press) == W22*roface
        // at every face, so the predictor hydrostatic force is zero -- including the surface band.
        // (The previous version built press from a per-cell analytic unit-gradient extrapolation,
        // whose face density disagreed with roface where the reinitialised phi is not exactly
        // unit-gradient -> a ~0.18 surface seed the multi-level coupling then amplified.)
        LEVEL_LOOP
        TILE_LOOP
        ILOOP
        JLOOP
        {
            const double dz  = p->amrex_geometry[p->level].CellSize(2);
            const double zlo = p->amrex_geometry[p->level].ProbLo(2);

            // --- absolute anchor for the tile-bottom cell (k=0 local) --------------------------
            // Integrate from the global column bottom with analytic unit-gradient phi (no off-tile
            // field reads -> order-independent). Exact in deep water; any small extrapolation error
            // near the surface only shifts the column's constant offset, which wpgrad never sees.
            k = 0;
            {
                const int    gk   = k + p->amr_tile_lo.z;    // GLOBAL k on this level
                const double zc   = p->pos_z();              // tile-bottom cell centre z
                const double phic = a->phi(i,j,k);
                const double zc0  = zlo + 0.5*dz;
                const double phi0 = phic + (zc - zc0);       // phi at global bottom, unit gradient
                double press = phi0*p->W1*fabs(p->W22) + p->I55;
                for(int m=0; m<gk; ++m)
                {
                    const double zf  = zlo + double(m+1)*dz;
                    const double phif= phic + (zc - zf);
                    const double H   = heaviside_ls(phif,psi);
                    const double rof = p->W1*H + p->W3*(1.0-H);
                    press += p->W22*dz*rof;                  // W22 < 0
                }
                a->press0(i,j,k) = a->press(i,j,k) = press;
            }

            // --- field-consistent upward integration (matches wpgrad's roface exactly) ---------
            for(k=0; k<KMAX_LOOP; ++k)
            {
                const double rof = dens.roface(p,a,0,0,1);   // face between cells k and k+1
                a->press0(i,j,k+1) = a->press(i,j,k+1)
                    = a->press(i,j,k) + p->W22*p->DZP[KP]*rof;
            }
        }

        // REEF_HYDRO_PROBE: measure the wpgrad-style vertical imbalance of the just-built press,
        // BEFORE and AFTER start4. If post-build ~0 but post-start4 != 0, the C-F interpolation in
        // start4 (FillPatchTwoLevels) is what breaks the hydrostatic balance, not the build.
        auto hydro_imb_report = [&](const char* tag)
        {
            if(!std::getenv("REEF_HYDRO_PROBE")) return;
            double worst=0.0; int wl=-1,wi[3]={-1,-1,-1}; double wpg=0,wrof=0;
            LEVEL_LOOP TILE_LOOP ILOOP JLOOP
            for(k=0;k<KMAX_LOOP;++k)
            {
                const double rof = dens.roface(p,a,0,0,1);
                const double dp  = a->press(i,j,k+1)-a->press(i,j,k);
                const double gg  = (a->grav_pot(i,j,k+1)-a->grav_pot(i,j,k))/p->DZP[KP];
                const double imb = std::fabs(-dp/(p->DZP[KP]*rof) + gg);
                if(imb>worst){ worst=imb; wl=p->level;
                    wi[0]=i+p->amr_tile_lo.x; wi[1]=j+p->amr_tile_lo.y; wi[2]=k+p->amr_tile_lo.z;
                    wpg=dp/(p->DZP[KP]*rof); wrof=rof; }
            }
            if(p->mpirank==0)
                std::cout<<"  [hydrobuild "<<tag<<"] worst |imb|="<<worst<<" at lev="<<wl
                         <<" ("<<wi[0]<<","<<wi[1]<<","<<wi[2]<<")  pgrad/rof="<<wpg
                         <<" roface="<<wrof<<std::endl;
        };
        hydro_imb_report("post-build");

        // REEF_COVERED_PRESS_RECON: match the per-step fix in pjm_corr::start -- overwrite covered
        // coarse press with the fine average so the iter-1 predictor already sees a horizontally
        // consistent covered column (the per-column hydrostatic offset otherwise seeds the geyser
        // from step 1). average_down writes only covered cells; uncovered columns untouched.
        if(std::getenv("REEF_COVERED_PRESS_RECON") && p->nlevs>1)
            for(int lev=p->nlevs-2; lev>=0; --lev)
            {
                amrex::average_down(a->press.GetMultiFab(lev+1),  a->press.GetMultiFab(lev),  0,1,p->ref_vec);
                amrex::average_down(a->press0.GetMultiFab(lev+1), a->press0.GetMultiFab(lev), 0,1,p->ref_vec);
            }

        // gcv 42 = transverse-linear C-F ghost (REEF_CF_GHOST_TRANSVERSE): reconstruct the coarse
        // press at the fine ghost's transverse position so the C-F ghost fill no longer breaks the
        // hydrostatic balance (post-start4 |imb| 0.22 -> ~0). Single-level falls back to gcv 40.
        const int gcv_h = (p->nlevs>1 && p->Y11==1) ? 42 : 40;
        pgc->start4(p,a->press,gcv_h,false);   // fill halos / C-F ghosts; NO average_down (keeps the
        pgc->start4(p,a->press0,gcv_h,false);  // coarse hydrostatic column self-consistent at surface)

        hydro_imb_report("post-start4");

        // REEF_HYDRO_INIT_PROBE: do coarse and fine hydrostatic press AGREE where they overlap?
        // Average the fine press down onto the covered coarse cells and compare with the directly
        // computed coarse press. The per-level dz quadrature of the SMOOTHED density (Hface over
        // psi) is dz-dependent across the surface, so a nonzero diff localised at the surface
        // (phi~0) confirms the levels carry inconsistent hydrostatic pressure there -- which the
        // C-F coupling then reconciles, breaking the fine surface balance.
        if(std::getenv("REEF_HYDRO_INIT_PROBE") && p->nlevs>1)
        {
            auto& cmf = a->press.GetMultiFab(0);
            amrex::MultiFab avgf(cmf.boxArray(), cmf.DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(avgf, cmf, 0, 0, 1, 0);
            amrex::average_down(a->press.GetMultiFab(1), avgf, 0, 1, p->ref_vec);

            auto& phimf = a->phi.GetMultiFab(0);
            double worst=0.0; int wi[3]={-1,-1,-1}; double w_pc=0,w_pa=0,w_phi=0;
            for(amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                const auto pc = cmf.const_array(mfi);
                const auto pa = avgf.const_array(mfi);
                const auto ph = phimf.const_array(mfi);
                for(int kk=bx.smallEnd(2); kk<=bx.bigEnd(2); ++kk)
                for(int jj=bx.smallEnd(1); jj<=bx.bigEnd(1); ++jj)
                for(int ii=bx.smallEnd(0); ii<=bx.bigEnd(0); ++ii)
                {
                    const double d = std::fabs(pc(ii,jj,kk) - pa(ii,jj,kk));
                    if(d > worst)
                    { worst=d; wi[0]=ii; wi[1]=jj; wi[2]=kk;
                      w_pc=pc(ii,jj,kk); w_pa=pa(ii,jj,kk); w_phi=ph(ii,jj,kk); }
                }
            }
            if(p->mpirank==0)
                std::cout<<"  [hydroinit] worst |press_coarse - avg(press_fine)|="<<worst
                         <<" at coarse ("<<wi[0]<<","<<wi[1]<<","<<wi[2]<<")"
                         <<"  press_c="<<w_pc<<" avgfine="<<w_pa<<" phi="<<w_phi
                         <<"  (nonzero @ surface => dz-inconsistent hydrostatic init)"<<std::endl;
        }
        #endif
    }

    if(p->I56==1)
    BASELOOP
    {
        if(a->phi(i,j,k)<0.0)
        a->press(i,j,k)=0.0;
    }
	
    BASELOOP
    a->press(i,j,k)+=p->I55;
	
	if(p->I12==2 && p->I30==0)
    GC4LOOP
	{
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];

        a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22) + p->I55;
	}
}
