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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"patchBC_interface.h"

void ioflow_f::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    p->gcin.resize_levels(nlevs);
    p->gcout.resize_levels(nlevs);
    int cs = 0;

    LEVEL_LOOP
    GC4LOOP
    {
        GCB4_TILE(n);

        i = p->gcb4[p->level][n].i;
        j = p->gcb4[p->level][n].j;
        k = p->gcb4[p->level][n].k;
        cs = p->gcb4[p->level][n].cs;

        if((p->gcb4[p->level][n].bc==1) && p->DF(i,j,k)>0)
        p->gcin[p->level].push_back({i, j, k, cs});
        else if((p->gcb4[p->level][n].bc==2) && p->DF(i,j,k)>0)
        p->gcout[p->level].push_back({i, j, k, cs});
    }
    GC_TILE_RESET;

    p->gcin_count=p->gcin[0].size();
    p->gcout_count=p->gcout[0].size();

    if(p->I10==1 && p->count==0)
    velini(p,a,pgc);


    // IO update

    GC4LOOP
    {
        GCB4_TILE(n);

        if(p->gcb4[p->level][n].bc==1 || p->gcb4[p->level][n].bc==6)
        {
            i = p->gcb4[p->level][n].i;
            j = p->gcb4[p->level][n].j;
            k = p->gcb4[p->level][n].k;

            if(p->DF(i,j,k)>0)
            {
                // inflow
                if(p->gcb4[p->level][n].cs==X_NEG)
                    p->IO[Im1JK] = 1;
                else if(p->gcb4[p->level][n].cs==X_POS)
                    p->IO[Ip1JK] = 1;
                else if(p->gcb4[p->level][n].cs==Y_NEG)
                    p->IO[IJm1K] = 1;
                else if(p->gcb4[p->level][n].cs==Y_POS)
                    p->IO[IJp1K] = 1;
                else if(p->gcb4[p->level][n].cs==Z_NEG)
                    p->IO[IJKm1] = 1;
                else if(p->gcb4[p->level][n].cs==Z_POS)
                    p->IO[IJKp1] = 1;
            }
        }

        if(p->gcb4[p->level][n].bc==2 || p->gcb4[p->level][n].bc==7)
        {
            i = p->gcb4[p->level][n].i;
            j = p->gcb4[p->level][n].j;
            k = p->gcb4[p->level][n].k;

            if(p->DF(i,j,k)>0)
            {
                // outflow
                if(p->gcb4[p->level][n].cs==X_NEG)
                    p->IO[Im1JK] = 2;
                else if(p->gcb4[p->level][n].cs==X_POS)
                    p->IO[Ip1JK] = 2;
                else if(p->gcb4[p->level][n].cs==Y_NEG)
                    p->IO[IJm1K] = 2;
                else if(p->gcb4[p->level][n].cs==Y_POS)
                    p->IO[IJp1K] = 2;
                else if(p->gcb4[p->level][n].cs==Z_NEG)
                    p->IO[IJKm1] = 2;
                else if(p->gcb4[p->level][n].cs==Z_POS)
                    p->IO[IJKp1] = 2;
            }
        }
    }
    GC_TILE_RESET;

    for(int qq=0;qq<pBC->obj_count;++qq)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    {
        if(pBC->patch[qq]->gcb[n][3]==1)
        p->IO[Im1JK] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==4)
        p->IO[Ip1JK] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==3)
        p->IO[IJm1K] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==2)
        p->IO[IJp1K] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==5)
        p->IO[IJKm1] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==6)
        p->IO[IJKp1] = 1;
    }
}

void ioflow_f::gcio_update_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    p->gcin.resize_levels(nlevs);
    p->gcout.resize_levels(nlevs);
    int cs = 0;

    LEVEL_LOOP
    GC4LOOP
    {
        GCB4_TILE(n);

        i = p->gcb4[p->level][n].i;
        j = p->gcb4[p->level][n].j;
        k = p->gcb4[p->level][n].k;
        cs = p->gcb4[p->level][n].cs;

        if((p->gcb4[p->level][n].bc==1) && p->DF(i,j,k)>0)
        p->gcin[p->level].push_back({i, j, k, cs});
        else if((p->gcb4[p->level][n].bc==2) && p->DF(i,j,k)>0)
        p->gcout[p->level].push_back({i, j, k, cs});
    }
    GC_TILE_RESET;

    p->gcin_count=p->gcin[0].size();
    p->gcout_count=p->gcout[0].size();

    // IO update
    MALOOP
    p->IO[IJK] = 0;

    GC4LOOP
    {
        GCB4_TILE(n);

        if(p->gcb4[p->level][n].bc==1 || p->gcb4[p->level][n].bc==6)
        {
            i = p->gcb4[p->level][n].i;
            j = p->gcb4[p->level][n].j;
            k = p->gcb4[p->level][n].k;

            if(p->DF(i,j,k)>0)
            {
                // inflow
                if(p->gcb4[p->level][n].cs==X_NEG)
                p->IO[Im1JK] = 1;
                else if(p->gcb4[p->level][n].cs==X_POS)
                p->IO[Ip1JK] = 1;
                else if(p->gcb4[p->level][n].cs==Y_NEG)
                p->IO[IJm1K] = 1;
                else if(p->gcb4[p->level][n].cs==Y_POS)
                p->IO[IJp1K] = 1;
                else if(p->gcb4[p->level][n].cs==Z_NEG)
                p->IO[IJKm1] = 1;
                else if(p->gcb4[p->level][n].cs==Z_POS)
                p->IO[IJKp1] = 1;
            }
        }

        if((p->gcb4[p->level][n].bc==2 || p->gcb4[p->level][n].bc==7))
        {
            i = p->gcb4[p->level][n].i;
            j = p->gcb4[p->level][n].j;
            k = p->gcb4[p->level][n].k;

            if(p->DF(i,j,k)>0)
            {
                // outflow
                if(p->gcb4[p->level][n].cs==X_NEG)
                p->IO[Im1JK] = 2;
                else if(p->gcb4[p->level][n].cs==X_POS)
                p->IO[Ip1JK] = 2;
                else if(p->gcb4[p->level][n].cs==Y_NEG)
                p->IO[IJm1K] = 2;
                else if(p->gcb4[p->level][n].cs==Y_POS)
                p->IO[IJp1K] = 2;
                else if(p->gcb4[p->level][n].cs==Z_NEG)
                p->IO[IJKm1] = 2;
                else if(p->gcb4[p->level][n].cs==Z_POS)
                p->IO[IJKp1] = 2;
            }
        }
    }
    GC_TILE_RESET;

    for(int qq=0;qq<pBC->obj_count;++qq)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    {
        if(pBC->patch[qq]->gcb[n][3]==1)
        p->IO[Im1JK] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==4)
        p->IO[Ip1JK] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==3)
        p->IO[IJm1K] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==2)
        p->IO[IJp1K] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==5)
        p->IO[IJKm1] = 1;
        else if(pBC->patch[qq]->gcb[n][3]==6)
        p->IO[IJKp1] = 1;
    }
}

void ioflow_f::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, reini *preini, ioflow *pflow)
{
    #if USE_AMREX
    walldin.resize(p->nlevs);
    #else
    walldin.resize(1);
    #endif
    // LEVEL_LOOP
    {
        walldin[p->level].resize(p->gcin.ssize(p->level));
        GCINLOOP
        {
            i=p->gcin[p->level][n].i;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            GCIN_TILE(n);

            walldin[p->level][n] = a->walld(i,j,k);
        }
        GC_TILE_RESET;
    }

    #if USE_AMREX
    walldout.resize(p->nlevs);
    #else
    walldout.resize(1);
    #endif
    // LEVEL_LOOP
    {
        walldout[p->level].resize(p->gcout.ssize(p->level));
        GCOUTLOOP
        {
            i=p->gcout[p->level][n].i;
            j=p->gcout[p->level][n].j;
            k=p->gcout[p->level][n].k;

            GCOUT_TILE(n);

            walldout[p->level][n] = a->walld(i,j,k);
        }
        GC_TILE_RESET;
    }
}

void ioflow_f::veltimesave(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    pvrans->veltimesave(p,a,pgc);
}

void ioflow_f::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc, vrans *pvrans)
{
    pvrans->sed_update(p,a,pgc);
}
