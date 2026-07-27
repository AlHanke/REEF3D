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
#include"fdm.h"

void ghostcell::set_DF(lexer *p, fdm *a)
{
    LOOP
    {
        int df = 1;

        if(p->solidread>0 && a->solid(i,j,k)<0.0)
            df = -1;
        else if(p->toporead>0 && a->topo(i,j,k)<0.0)
            df = -1;
        else if(p->X10>0 && a->fb(i,j,k)<0.0)
            df = -1;

        p->DF(i,j,k) = df;
    }

    #if USE_AMREX
    p->DF.fillBoundary();
    #else
    flagx(p,p->DF);
    #endif
}

void ghostcell::gcdf_update(lexer *p, fdm *a)
{
    // -----------------------------------------------------------
    // FLAG
    if(p->G5==1)
    {
        double psi;
        p->level = 0;
        TILE_LOOP
        IJKLOOP
        PBASECHECK
        {
            if (p->j_dir==0)
                psi = -p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
            else if (p->j_dir==1)
                psi = -p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

            if((a->solid(i,j,k)>=psi || p->solidread==0) && (a->topo(i,j,k)>=psi || p->toporead==0))
                p->flag4[IJK]=WATER_FLAG;

            if((a->solid(i,j,k)<psi && p->solidread==1) || (a->topo(i,j,k)<psi && p->toporead==1))
                p->flag4[IJK]=-10;
        }

        #if USE_AMREX
        p->flag4.fillBoundary();
        #else
        flagx(p,p->flag4);
        #endif

        p->level = 0;
        TILE_LOOP
        IJKLOOP
        PBASECHECK
        {
            p->flag1[IJK]=p->flag4[IJK];
            p->flag2[IJK]=p->flag4[IJK];
            p->flag3[IJK]=p->flag4[IJK];

            if(p->flag4[IJK]>0 && p->flag4[Ip1JK]<0)
                p->flag1[IJK]=-10;

            if(p->flag4[IJK]>0 && p->flag4[IJp1K]<0)
                p->flag2[IJK]=-10;

            if(p->flag4[IJK]>0 && p->flag4[IJKp1]<0)
                p->flag3[IJK]=-10;
        }

        gcb_velflagio(p);

        #if USE_AMREX
        p->flag1.fillBoundary();
        p->flag2.fillBoundary();
        p->flag3.fillBoundary();
        #else
        flagx(p,p->flag1);
        flagx(p,p->flag2);
        flagx(p,p->flag3);
        #endif

        #if USE_AMREX
        p->flag1.fillHigherLevels();
        p->flag2.fillHigherLevels();
        p->flag3.fillHigherLevels();
        p->flag4.fillHigherLevels();
        #endif
    }

    // -----------------------------------------------------------
    // DF
    set_DF(p,a);

    // -----------------------------------------------------------
    // count gcdf entries
    count=0;

    // gcdf count
    p->level = 0;
    TILE_LOOP
    IJKLOOP
    PBASECHECK
    if(p->DF(i,j,k)>0)
    {
        if(p->DF(i-1,j,k)<0)
            ++count;

        if(p->DF(i+1,j,k)<0)
            ++count;

        if(p->DF(i,j-1,k)<0)
            ++count;

        if(p->DF(i,j+1,k)<0)
            ++count;

        if(p->DF(i,j,k-1)<0)
            ++count;

        if(p->DF(i,j,k+1)<0)
            ++count;
    }

    if(p->gcdf4_count!=count)
    {
        p->gcdf4.resize(count);

        p->gcdf4_count=count;
    }

    // assign gcdf entries
    //
    // (i,j,k) recorded here are TILE-LOCAL — they are offsets from the producing
    // tile's amr_tile_lo. Consumers run outside TILE_LOOP, where the default
    // whole-box context is active, so each entry also records the dense id
    // (column 6) of the tile it was written under. Without it the indices resolve
    // against the wrong origin as soon as tiling produces more than one tile per
    // box. TILE_CTX_DEFAULT here would mean the entry was produced outside a tile
    // loop; it is not, but the id records that faithfully either way.
    count=0;

    p->level = 0;
    TILE_LOOP
    {
        #if USE_AMREX
        const int tilehandle = p->tile_ctx_id();
        #else
        const int tilehandle = 0;
        #endif

        IJKLOOP
        PBASECHECK
        if(p->DF(i,j,k)>0)
        {
            if(p->DF(i-1,j,k)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=1;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }

            if(p->DF(i+1,j,k)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=4;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }

            if(p->DF(i,j-1,k)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=3;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }

            if(p->DF(i,j+1,k)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=2;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }

            if(p->DF(i,j,k-1)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=5;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }

            if(p->DF(i,j,k+1)<0)
            {
                p->gcdf4[count][0]=i;
                p->gcdf4[count][1]=j;
                p->gcdf4[count][2]=k;
                p->gcdf4[count][3]=6;
                p->gcdf4[count][5]=tilehandle;
                ++count;
            }
        }
    }

    fieldint4 cval(p);

    count=0;

    p->level = 0;
    TILE_LOOP
    IJKLOOP
    PBASECHECK
    {
        cval(i,j,k)=count;

        ++count;
    }

    // gcb4 entries carry the tile they were generated under, so this resolves
    // the same way the gcdf4 loop below does.
    GC4LOOP
    {
        GCB4_TILE(n);

        i=p->gcb4[p->level][n].i;
        j=p->gcb4[p->level][n].j;
        k=p->gcb4[p->level][n].k;
        p->gcb4[p->level][n].row=cval(i,j,k);
    }
    GC_TILE_RESET;


    GCDF4LOOP
    {
        GCDF4_TILE(n);
        i=p->gcdf4[n][0];
        j=p->gcdf4[n][1];
        k=p->gcdf4[n][2];
        p->gcdf4[n][4]=cval(i,j,k);
    }
    GC_TILE_RESET;

    // flagsf1/2/3's replacements are the staggered-face versions of DF (same
    // convention as flag1-3 vs flag4): open unless the face's neighbor cell is
    // blocked. Since that's a pure function of DF -- which already carries a
    // p->margin-deep synced halo -- they're computed on the fly here instead of
    // as separate arrays; gcdf_update_impl never reaches more than 2 cells from DF.
    auto flagsf1 = [p](int i, int j, int k) -> int
    {
        int v = p->DF(i,j,k);
        return (v>0 && p->DF(i+1,j,k)<0) ? -1 : v;
    };

    auto flagsf2 = [p](int i, int j, int k) -> int
    {
        int v = p->DF(i,j,k);
        return (v>0 && p->DF(i,j+1,k)<0) ? -1 : v;
    };

    auto flagsf3 = [p](int i, int j, int k) -> int
    {
        int v = p->DF(i,j,k);
        return (v>0 && p->DF(i,j,k+1)<0) ? -1 : v;
    };

    gcdf_update_impl(p, flagsf1, p->gcdf1, p->gcdf1_count);
    gcdf_update_impl(p, flagsf2, p->gcdf2, p->gcdf2_count);
    gcdf_update_impl(p, flagsf3, p->gcdf3, p->gcdf3_count);
}

template<typename FlagT, typename GcdfT>
void ghostcell::gcdf_update_impl(lexer *p, FlagT &flagsf, GcdfT &gcdf, int &gcdf_count)
{
    // -----------------------
    // flagsf

    count = 0;

    p->level = 0;
    TILE_LOOP
    IJKLOOP
    PBASECHECK
    if(flagsf(i,j,k)>0)
    {
        if(flagsf(i-1,j,k)<0)
            ++count;

        if(flagsf(i+1,j,k)<0)
            ++count;

        if(flagsf(i,j-1,k)<0)
            ++count;

        if(flagsf(i,j+1,k)<0)
            ++count;

        if(flagsf(i,j,k-1)<0)
            ++count;

        if(flagsf(i,j,k+1)<0)
            ++count;
    }

    if(gcdf_count!=count)
    {
        gcdf.assign(count, {});

        gcdf_count=count;
    }

    // assign gcdf entries (tile-local i,j,k + tile handle, see gcdf4 above)
    count=0;

    p->level = 0;
    TILE_LOOP
    {
        #if USE_AMREX
        const int tilehandle = p->tile_ctx_id();
        #endif

        IJKLOOP
        PBASECHECK
        if(flagsf(i,j,k)>0)
        {
            if(flagsf(i-1,j,k)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=1;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }

            if(flagsf(i+1,j,k)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=4;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }

            if(flagsf(i,j-1,k)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=3;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }

            if(flagsf(i,j+1,k)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=2;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }

            if(flagsf(i,j,k-1)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=5;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }

            if(flagsf(i,j,k+1)<0)
            {
                gcdf[count].i=i;
                gcdf[count].j=j;
                gcdf[count].k=k;
                gcdf[count].cs=6;
                #if USE_AMREX
                gcdf[count].ctx_id=tilehandle;
                #endif
                ++count;
            }
        }
    }
}
