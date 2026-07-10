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
    TILE_LOOP
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

    p->level = 0;
    TILE_LOOP
    MALOOP
    {
        p->flag1(i,j,k)=p->flag4(i,j,k);
        p->flag2(i,j,k)=p->flag4(i,j,k);
        p->flag3(i,j,k)=p->flag4(i,j,k);
    }

    GC4LOOP
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        if(p->gcb4[n][3]==4 && (p->periodic1!=1 || i+p->origin_i<p->gknox-1))
            p->flag1[IJK]=OBJ_FLAG;
    }

    GC4LOOP
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        if(p->gcb4[n][3]==2 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
            p->flag2[IJK]=OBJ_FLAG;
    }

    GC4LOOP
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        if(p->gcb4[n][3]==6 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
            p->flag3[IJK]=OBJ_FLAG;
    }

    #if USE_AMREX
    p->flag1.fillHigherLevels();
    p->flag2.fillHigherLevels();
    p->flag3.fillHigherLevels();
    p->flag4.fillHigherLevels();
    #else
    flagx(p,p->flag1);
    flagx(p,p->flag2);
    flagx(p,p->flag3);
    #endif

    #if USE_AMREX
    // Bottom solid-wall ghosts -> OBJ_FLAG. The -1->OBJ conversion above runs before
    // fillBoundary(), which overwrites the domain-boundary ghost ring, leaving the bottom wall
    // ghost as an AIR-like -1 on the coarse level -- and fillHigherLevels leaves the FINE-level
    // exterior ghost at 0 (uninitialised). poisson_pcorr's wall fold needs flag4<0, so both an
    // AIR-like -1 AND a 0 ghost skip the fold, leaving a spurious floor coupling in the matrix
    // (projcheck operator residual O(1e-3) at fine floor cells). Re-tag ONLY the z-low (bottom)
    // ghost ring -- the single true solid wall of this domain. Condition catches BOTH f==0 (fine)
    // and f<0 air-like (coarse); the lateral/top boundaries are symmetry planes (kept), the free
    // surface is internal (kept Dirichlet), INFLOW/OUTFLOW and already-solid are preserved.
    LEVEL_LOOP
    {
        const auto dom = p->amrex_geometry[p->level].Domain();
        TILE_LOOP
        for (i = -margin; i <= (p->amr_tile_hi.x - p->amr_tile_lo.x)+margin; ++i)
        for (j = -margin; j <= (p->amr_tile_hi.y - p->amr_tile_lo.y)+margin; ++j)
        for (k = -margin; k <= (p->amr_tile_hi.z - p->amr_tile_lo.z)+margin; ++k)
        {
            const int gk = k + p->amr_tile_lo.z;
            const int f  = p->flag4(i,j,k);
            if(gk<dom.smallEnd(2) && f!=INFLOW_FLAG && f!=OUTFLOW_FLAG && f>SOLID_FLAG && f<WATER_FLAG)
                p->flag4(i,j,k)=OBJ_FLAG;
        }
    }
    #endif
}
