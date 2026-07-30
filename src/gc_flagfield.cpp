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

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::flagfield(lexer *p)
{
    LEVEL_LOOP
    TILE_LOOP
    MALOOP
    {
        if(p->flag4(i,j,k)==1)
        {
            p->flag4(i,j,k)=WATER_FLAG;
        }
        else if(p->flag4(i,j,k)==-1)
        {
            p->flag4(i,j,k)=OBJ_FLAG;
        }
    }

    #if USE_AMREX
    p->flag4.fillBoundary();
    #else
    flagx(p,p->flag4);
    #endif

    p->level = 0;
    if(p->Y60==1)
    // TILE_LOOP
    IJKLOOP
    PCHECK
    {
        if(p->i_dir==1)
        if(p->flag4[Im1JK]<0
        && p->flag4[Ip1JK]<0)
            p->flag4[IJK]=OBJ_FLAG;

        if(p->j_dir==1)
        if(p->flag4[IJm1K]<0
        && p->flag4[IJp1K]<0)
            p->flag4[IJK]=OBJ_FLAG;

        if(p->k_dir==1)
        if(p->flag4[IJKm1]<0
        && p->flag4[IJKp1]<0)
            p->flag4[IJK]=OBJ_FLAG;
    }

    #if USE_AMREX
    p->flag4.fillHigherLevels();
    #endif

    #if USE_AMREX
    // Domain-boundary ghost flag4 on every level -> OBJ_FLAG at CLOSED walls. The -1->OBJ
    // conversion above runs before fillBoundary(), which overwrites the domain-boundary ghost
    // ring, leaving the ghost AIR-like (-1) on the coarse level -- and fillHigherLevels leaves the
    // FINE-level exterior ghost at 0 (uninitialised). poisson_pcorr's wall-Neumann fold (solidnb)
    // only fires for a solid flag4<0, so a 0 (fine) OR air-like -1 (coarse) ghost at a closed
    // boundary SKIPS the fold, leaving a dangling coupling HYPRE drops -> a non-conservative row ->
    // projection blow-up. Static built all its walls at init (via the full flag machinery) so was
    // immune; an adaptive regrid only re-derives the fine level through fillHigherLevels, so its
    // lateral/top walls arrive as 0 -- the "AMR stop after 1 iteration".
    //
    // Re-tag every CLOSED, non-periodic domain-boundary ghost. A solid WALL (bcside==21) and a
    // SYMMETRY plane (bcside==3) both take homogeneous Neumann on the pressure matrix, so both map
    // to OBJ_FLAG; INFLOW (1,6) / OUTFLOW (2,7) are excluded (they carry INFLOW/OUTFLOW_FLAG and a
    // different BC). Velocity BCs are unaffected -- they come from the bcside-driven BCRec fill in
    // start1/2/3, not flag4. Face->bcside map (lexer, per grid_amrex bc_type): x-/x+ =1/4,
    // y-/y+ =3/2, z-/z+ =5/6. Already-solid, water/interior, INFLOW and OUTFLOW cells are preserved.
    // CRITICAL level-0 caveat: the z-low (bottom) wall is re-tagged on EVERY level (it is a true
    // solid wall and the coarse bottom arrives as air-like -1). The lateral/top walls, however, are
    // re-tagged on the FINE levels ONLY. On the coarse level the free-surface Dirichlet anchor is
    // COVERED away wherever a finer level refines the interface (identity rows), so the lateral/top
    // walls staying air/Dirichlet are what keep the coarse pressure operator non-singular; forcing
    // them Neumann there makes the (surface-covered) coarse solve singular -> the solver diverges
    // (pres ~1e10) and the run overflows. The fine level carries the free surface itself, so its
    // walls are safely homogeneous Neumann.
    auto closed_bc = [](int bc){ return bc!=1 && bc!=2 && bc!=6 && bc!=7; }; // wall(21)/symmetry(3): fold
    const bool xlo_c = closed_bc(p->bcside1) && p->periodic1==0;
    const bool xhi_c = closed_bc(p->bcside4) && p->periodic1==0;
    const bool ylo_c = closed_bc(p->bcside3) && p->periodic2==0;
    const bool yhi_c = closed_bc(p->bcside2) && p->periodic2==0;
    const bool zlo_c = closed_bc(p->bcside5) && p->periodic3==0;
    const bool zhi_c = closed_bc(p->bcside6) && p->periodic3==0;

    LEVEL_LOOP
    {
        const bool fine = (p->level > 0);
        LEVEL_DOMAIN_DECL(dom)
        TILE_LOOP
        for (i = -margin; i <= (p->amr_tile_hi.x - p->amr_tile_lo.x)+margin; ++i)
        for (j = -margin; j <= (p->amr_tile_hi.y - p->amr_tile_lo.y)+margin; ++j)
        for (k = -margin; k <= (p->amr_tile_hi.z - p->amr_tile_lo.z)+margin; ++k)
        {
            const int f = p->flag4(i,j,k);
            if(f==INFLOW_FLAG || f==OUTFLOW_FLAG || f<=SOLID_FLAG || f>=WATER_FLAG) continue;

            const int gi = GLOBAL_I;
            const int gj = GLOBAL_J;
            const int gk = GLOBAL_K;

            const bool wall =
                   (gk<DOMAIN_LO(dom,2) && zlo_c)                                   // z-low: ALL levels
                || (fine && ( (gi<DOMAIN_LO(dom,0) && xlo_c) || (gi>DOMAIN_HI(dom,0) && xhi_c)
                           || (gj<DOMAIN_LO(dom,1) && ylo_c) || (gj>DOMAIN_HI(dom,1) && yhi_c)
                           || (gk>DOMAIN_HI(dom,2) && zhi_c) ));                     // lateral/top: FINE only
            if(wall)
                p->flag4(i,j,k)=OBJ_FLAG;
        }
    }
    p->flag4.fillBoundary();
    #else
    flagx(p,p->flag4);
    #endif

    // gcb4 is derived from flag4, not read from the grid file. It has to be
    // rebuilt here rather than once at setup: the loops below consume it, and
    // regrid re-chops the decomposition the entries' indices are relative to.
    // Placed after the OBJ_FLAG conversion and FillBoundary above (so the
    // domain-exterior ghost ring reads OBJ_FLAG) and before the tagging below.
    gcb4_generate(p);

    LEVEL_LOOP
    TILE_LOOP
    MALOOP
    {
        p->flag1(i,j,k)=p->flag4(i,j,k);
        p->flag2(i,j,k)=p->flag4(i,j,k);
        p->flag3(i,j,k)=p->flag4(i,j,k);
    }

    // Domain-high guard for the periodic directions, in the LEVEL's index space —
    // GLOBAL_I/J/K and DOMAIN_HI from global_index.h, which documents why the legacy
    // i+origin_i < gknox-1 form is wrong in both operands under AMReX. DOMAIN_HI is
    // gknox-1 at level 0, so level-0 behaviour is unchanged.
    //
    // Belt-and-braces on the AMReX path: FillBoundary now wraps periodic directions,
    // so a periodic domain face reads the interior value through its ghost, no entry
    // is emitted there by gcb4_generate, and this guard never fires. If it ever does,
    // the finding is that the periodic ghost was not filled.
    LEVEL_LOOP
    {
        // Read before any GCB4_TILE, i.e. at the loop's own level. Safe against the
        // set_tile_ctx inside the body, which restores the entry's level — equal to
        // the bucket the entry lives in, hence to this one.
        LEVEL_DOMAIN_DECL(dom)

        GC4LOOP
        {
            GCB4_TILE(n);

            i=p->gcb4[p->level][n].i;
            j=p->gcb4[p->level][n].j;
            k=p->gcb4[p->level][n].k;

            if(p->gcb4[p->level][n].cs==X_POS && (p->periodic1!=1 || GLOBAL_I<DOMAIN_HI(dom,0)))
                p->flag1(i,j,k)=OBJ_FLAG;

            if(p->gcb4[p->level][n].cs==Y_POS && (p->periodic2!=1 || GLOBAL_J<DOMAIN_HI(dom,1)))
                p->flag2(i,j,k)=OBJ_FLAG;

            if(p->gcb4[p->level][n].cs==Z_POS && (p->periodic3!=1 || GLOBAL_K<DOMAIN_HI(dom,2)))
                p->flag3(i,j,k)=OBJ_FLAG;
        }
    }
    GC_TILE_RESET;

    // Propagate the tagging above into the same-level ghost ring: GC4LOOP writes
    // OBJ_FLAG only at gcb4 cells, which are valid cells, so a neighbouring box's
    // ghost copy of a tagged cell would otherwise stay pre-tag and any stencil
    // reading flag1/2/3 across a box seam would see the stale value.
    //
    // AMReX-only: ArrayWrapper3D::fillBoundary is compiled only under USE_AMREX.
    // The legacy counterpart would be flagx(p,p->flag1) etc. — deliberately NOT
    // added, since flagfield has never synced flag1/2/3 here in that build and this
    // is not the place to change its behaviour untested.
    #if USE_AMREX
    p->flag1.fillBoundary();
    p->flag2.fillBoundary();
    p->flag3.fillBoundary();
    #else
    flagx(p,p->flag1);
    flagx(p,p->flag2);
    flagx(p,p->flag3);
    #endif
}
