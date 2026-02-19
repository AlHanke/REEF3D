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
    p->define_inflow_outflow_ba();
    #else
    int count1,count2;

    count1=0;
    count2=0;
    GC4LOOP
    {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];

        if((p->gcb4[n][4]==1 || p->gcb4[n][4]==6) && p->DF[IJK]>0)
        ++count1;

        if((p->gcb4[n][4]==2 || p->gcb4[n][4]==7) && p->DF[IJK]>0)
        ++count2;
    }

    p->Iresize(p->gcin,p->gcin_count, count1, 4, 4);
    p->Iresize(p->gcout,p->gcout_count, count2, 4, 4);

    count1=0;
    count2=0;
    GC4LOOP
    {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];

        if(p->gcb4[n][4]==1 && p->DF[IJK]>0)
        {
            p->gcin[count1][0]=p->gcb4[n][0];
            p->gcin[count1][1]=p->gcb4[n][1];
            p->gcin[count1][2]=p->gcb4[n][2];
            p->gcin[count1][3]=p->gcb4[n][3];
            ++count1;
        }
        else if(p->gcb4[n][4]==2 && p->DF[IJK]>0)
        {
            p->gcout[count2][0]=p->gcb4[n][0];
            p->gcout[count2][1]=p->gcb4[n][1];
            p->gcout[count2][2]=p->gcb4[n][2];
            p->gcout[count2][3]=p->gcb4[n][3];
            ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
    #endif

    if(p->I10==1 && p->count==0)
    velini(p,a,pgc);


    // IO update

    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
            i = p->gcb4[n][0];
            j = p->gcb4[n][1];
            k = p->gcb4[n][2];

            if(p->DF[IJK]>0)
            {
                // inflow
                if(p->gcb4[n][3]==1)
                    p->IO[Im1JK] = 1;
                else if(p->gcb4[n][3]==4)
                    p->IO[Ip1JK] = 1;
                else if(p->gcb4[n][3]==3)
                    p->IO[IJm1K] = 1;
                else if(p->gcb4[n][3]==2)
                    p->IO[IJp1K] = 1;
                else if(p->gcb4[n][3]==5)
                    p->IO[IJKm1] = 1;
                else if(p->gcb4[n][3]==6)
                    p->IO[IJKp1] = 1;
            }
        }

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7)
        {
            i = p->gcb4[n][0];
            j = p->gcb4[n][1];
            k = p->gcb4[n][2];

            if(p->DF[IJK]>0)
            {
                // outflow
                if(p->gcb4[n][3]==1)
                    p->IO[Im1JK] = 2;
                else if(p->gcb4[n][3]==4)
                    p->IO[Ip1JK] = 2;
                else if(p->gcb4[n][3]==3)
                    p->IO[IJm1K] = 2;
                else if(p->gcb4[n][3]==2)
                    p->IO[IJp1K] = 2;
                else if(p->gcb4[n][3]==5)
                    p->IO[IJKm1] = 2;
                else if(p->gcb4[n][3]==6)
                    p->IO[IJKp1] = 2;
            }
        }
    }

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
    int count1,count2;

    count1=0;
    count2=0;
    GC4LOOP
    {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];

        if((p->gcb4[n][4]==1 || p->gcb4[n][4]==6) && p->DF[IJK]>0 && p->wet[IJ]==1)
        ++count1;

        if((p->gcb4[n][4]==2 || p->gcb4[n][4]==7) && p->DF[IJK]>0)
        ++count2;
    }

    p->Iresize(p->gcin,p->gcin_count, count1, 4, 4);
    p->Iresize(p->gcout,p->gcout_count, count2, 4, 4);

    count1=0;
    count2=0;
    GC4LOOP
    {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];

        if((p->gcb4[n][4]==1 || p->gcb4[n][4]==6) && p->DF[IJK]>0 && p->wet[IJ]==1)
        {
            p->gcin[count1][0]=p->gcb4[n][0];
            p->gcin[count1][1]=p->gcb4[n][1];
            p->gcin[count1][2]=p->gcb4[n][2];
            p->gcin[count1][3]=p->gcb4[n][3];
            ++count1;
        }

        if((p->gcb4[n][4]==2 || p->gcb4[n][4]==7) && p->DF[IJK]>0)
        {
            p->gcout[count2][0]=p->gcb4[n][0];
            p->gcout[count2][1]=p->gcb4[n][1];
            p->gcout[count2][2]=p->gcb4[n][2];
            p->gcout[count2][3]=p->gcb4[n][3];
            ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;

    // IO update
    MALOOP
    p->IO[IJK] = 0;

    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
            i = p->gcb4[n][0];
            j = p->gcb4[n][1];
            k = p->gcb4[n][2];

            if(p->DF[IJK]>0)
            {
                // inflow
                if(p->gcb4[n][3]==1)
                p->IO[Im1JK] = 1;
                else if(p->gcb4[n][3]==4)
                p->IO[Ip1JK] = 1;
                else if(p->gcb4[n][3]==3)
                p->IO[IJm1K] = 1;
                else if(p->gcb4[n][3]==2)
                p->IO[IJp1K] = 1;
                else if(p->gcb4[n][3]==5)
                p->IO[IJKm1] = 1;
                else if(p->gcb4[n][3]==6)
                p->IO[IJKp1] = 1;
            }
        }

        if((p->gcb4[n][4]==2 || p->gcb4[n][4]==7))
        {
            i = p->gcb4[n][0];
            j = p->gcb4[n][1];
            k = p->gcb4[n][2];

            if(p->DF[IJK]>0)
            {
                // outflow
                if(p->gcb4[n][3]==1)
                p->IO[Im1JK] = 2;
                else if(p->gcb4[n][3]==4)
                p->IO[Ip1JK] = 2;
                else if(p->gcb4[n][3]==3)
                p->IO[IJm1K] = 2;
                else if(p->gcb4[n][3]==2)
                p->IO[IJp1K] = 2;
                else if(p->gcb4[n][3]==5)
                p->IO[IJKm1] = 2;
                else if(p->gcb4[n][3]==6)
                p->IO[IJKp1] = 2;
            }
        }
    }

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
    LEVEL_LOOP
    {
        #if USE_AMREX
        walldin[p->level].resize(p->inflow_ijk[p->level].size());
        for(n=0; n<p->inflow_ijk[p->level].size(); n++)
        {
            auto iv = p->inflow_ijk[p->level][n];
            i=iv[0];
            j=iv[1];
            k=iv[2];
        #else
        walldin[p->level].resize(p->gcin_count);
        for(n=0;n<p->gcin_count;++n)
        {
            i=p->gcin[n][0];
            j=p->gcin[n][1];
            k=p->gcin[n][2];
        #endif

            walldin[p->level][n] = a->walld(i,j,k);
        }
    }

    #if USE_AMREX
    walldout.resize(p->nlevs);
    #else
    walldout.resize(1);
    #endif
    LEVEL_LOOP
    {
        #if USE_AMREX
        walldout[p->level].resize(p->outflow_ijk[p->level].size());
        for(n=0; n<p->outflow_ijk[p->level].size(); n++)
        {
            auto iv = p->outflow_ijk[p->level][n];
            i=iv[0];
            j=iv[1];
            k=iv[2];
        #else
        walldout[p->level].resize(p->gcout_count);
        for(n=0;n<p->gcout_count;++n)
        {
            i=p->gcout[n][0];
            j=p->gcout[n][1];
            k=p->gcout[n][2];
        #endif

            walldout[p->level][n] = a->walld(i,j,k);
        }
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
