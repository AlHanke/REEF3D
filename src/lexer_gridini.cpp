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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#if USE_AMREX
#include "fieldint_amrex.h"
#include "reini.h"
#include "6DOF.h"
#include "fdm.h"
#endif

void lexer::gridini(ghostcell *pgc)
{
    #if USE_AMREX
    setup_amrex_geometry(this,pgc);
    #endif

    grid::gridspacing(this, pgc);
    #if USE_AMREX
    if(nlevs > 1)
    {
        if(G2==1)
        {
            std::cerr << "Error: G2=1 is not supported for nlevs > 1" << std::endl;
            pgc->final(true);
        }

        update_cell_coordinates();
        update_cell_spacing();
    }
    #endif
}

void lexer::flagini()
{
    control_calc();
	gridsize();

    lexer *p = this;

    flag4.resize(-1);

    p->level = 0;
    // TILE_LOOP
    IJKLOOP
    flag4(i,j,k) = flag4_grid[IJK];

    flag4_grid.reset(); // transported into flag4; not needed afterwards

    flag1.resize(OBJ_FLAG);
    flag2.resize(OBJ_FLAG);
    flag3.resize(OBJ_FLAG);
    flag5.resize();

    // boundary conditions
    IO.resize(0);
    IOSL.resize(0);
    DF.resize(1);
#if USE_AMREX
    m_df123 = make_imf(this, 3, &m_df123);
#endif
    // DF1/DF2/DF3: view mode under USE_AMREX (resize() is a no-op; data lives in m_df123)
    DF1.resize();
    DF2.resize();
    DF3.resize();

    // gcdf
    gcdf1_count=gcdf2_count=gcdf3_count=gcdf4_count=1;

    // gcsldf
    gcsldfeta4_count=1;

    gcsldfeta4.resize(gcsldfeta4_count);

    sliceflagini();
}

int lexer::conv(double a)
{
	int b,c;
	double d,diff;

	c= int( a);
	d=double(c);
	diff=a-d;

	b=c;

	if(diff>0.5)
	b=c+1;

	if(diff<=-0.5)
	b=c-1;

	return b;
}

void lexer::regrid(fdm* a, reini* preini, sixdof* p6dof, ghostcell* pgc, ioflow* pflow)
{
    #if USE_AMREX
    changed = false;
    // DIAGNOSTIC (Step 1 of plan fluffy-oasis): gate regrid to every 10 steps
    // to isolate per-step BoxArray churn + pc_interp dissipation. Revert when done.
    if (count % 10000000 == 0 && count > 0)
    {
        // Localisation probe: max|field| over valid cells per level at each sub-step of
        // regrid, so a blow-up can be traced to the exact stage that injects it (rebuild
        // vs reinit vs ghost restore). norm0 does its own MPI reduce. REEF_REGRID_PROBE.
        auto probe = [&](const char* tag)
        {
            if (!std::getenv("REEF_REGRID_PROBE")) return;
            double mu=0,mv=0,mw=0,mp=0,mphi=0;
            for (int lev=0; lev<nlevs; ++lev)
            {
                mu   = std::max(mu,   a->u.GetMultiFab(lev).norm0());
                mv   = std::max(mv,   a->v.GetMultiFab(lev).norm0());
                mw   = std::max(mw,   a->w.GetMultiFab(lev).norm0());
                mp   = std::max(mp,   a->press.GetMultiFab(lev).norm0());
                mphi = std::max(mphi, a->phi.GetMultiFab(lev).norm0());
            }
            if (mpirank==0)
                std::cout << "  [regridprobe] count=" << count << " nlevs=" << nlevs
                          << " " << tag << "  |u|=" << mu << " |v|=" << mv << " |w|=" << mw
                          << " |press|=" << mp << " |phi|=" << mphi << std::endl;
        };

        probe("entry");
        grid_amrex::regrid_amrex_box_array_and_distribution_mapping(this, a); // Bug with higher levels of static refinement
        grid_amrex::update_cell_coordinates();
        grid_amrex::update_cell_spacing();
        grid_amrex::update_registered_weno(nlevs);

        // Re-establish the flag semantics for the re-chopped grid. The field rebuild inside
        // regrid_amrex_box_array_and_distribution_mapping reallocates the registered flag
        // MultiFabs (flag1-4): valid cells are restored from the old data, but the ghost ring
        // is left at setVal(0) (an INVALID flag) and flag1/2/3 are not re-derived. gridini runs
        // flagini()+flagfield() once at setup and the grid never changed thereafter; regrid must
        // repeat flagfield() so the WATER/OBJ conversion, wall tagging (flag1/2/3 via gcb4),
        // fillBoundary and the C-F fillHigherLevels() are rebuilt. Without it the projection has
        // no wall BCs on the new grid -> it cannot balance the gravity predictor -> velocity
        // leaks and grows (post-mom |u| jumps from 0 to O(1) at rest).
        pgc->flagfield(this);
        gridhelper();
        // rebuild gcio
        probe("post-rebuild");

        // Restore phi's ghost band BEFORE reinit reads it. The field rebuild inside
        // regrid_amrex_box_array_and_distribution_mapping drops phi's physical/coarse-fine
        // ghosts (FillBoundary fills same-level/periodic only), and reini_RK3::start's first
        // prdisc->start reads phi(i+-1) with no preceding fill when count>0 -- so stale ghost
        // junk would be swept into valid phi by the redistance stencil. gcval_phi = 50+F50
        // matches reini_RK3's level-set boundary variant (F50==1->51 ... 4->54). At count==0
        // reinit self-fills phi's ghost (gcval_iniphi), so this is a no-op there.
        pgc->start4(this,a->phi,50+F50);
        probe("post-phi-ghost");

        preini->start(a,this,a->phi,pgc,pflow);
        probe("post-reinit");
        lexer* p = this;
        int counter = 0;
        PLAINLOOP
        {
            counter++;
        }
        veclength += counter - cellnum;
        cellnum = counter;
        a->rhsvec.resize(veclength);
        a->M.resize(veclength);
        LEVEL_LOOP
        TILE_LOOP
        MALOOP
        {
            a->grav_pot(i,j,k) = p->W20*p->pos_x() + p->W21*p->pos_y() + p->W22*p->pos_z();
        }

        // Restore the ghost band dropped by the field rebuild inside
        // regrid_amrex_box_array_and_distribution_mapping. fill_registered_mf_level
        // reallocates each registered MultiFab and reloads only the valid region
        // (InterpFromCoarseLevel / ParallelCopy fill zero ghost cells; the trailing
        // FillBoundary fills same-level and periodic ghosts only). Physical-boundary
        // and coarse-fine ghost cells are owned by ghostcell::start1-4 (FillDomainBoundary),
        // which regrid never calls -- so without this the freshly-allocated ghost band
        // is read as garbage by the next predictor (exceeding velocities) or, if a
        // rebalance follows immediately, by its grad(press)/roface stencil (HYPRE Inf/NaN).
        // phi is handled separately (pre-reinit fill above + preini->start's own re-sync).
        pgc->start1(p,a->u,10);
        pgc->start2(p,a->v,11);
        pgc->start3(p,a->w,12);
        pgc->start4(p,a->press,40,false);
        // pgc->gcdf_update(p,a);
        probe("post-ghost-restore (exit)");
    }
    #endif
}