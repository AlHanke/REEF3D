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
#include"field.h"
#if USE_AMREX
#include "field_amrex.h"
#endif

void ghostcell::start1(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    #if USE_AMREX
    starttime=timer();
    f.FillDomainBoundary(gcv);
    endtime=timer();
    p->xtime+=endtime-starttime;
    p->gctime+=endtime-starttime;
    #else
    if(do_comms)
    {
        starttime=timer();
        gcparax(p,f,1);
        gcparacox(p,f);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    // solid ghostcells
    starttime=timer();
    QQGC1LOOP
    if((p->gcb1[p->level][qq].cs!=2 && p->gcb1[p->level][qq].cs!=3) || p->j_dir==1)
        gcdistro1(f, p->gcb1[p->level][qq].i, p->gcb1[p->level][qq].j, p->gcb1[p->level][qq].k, p->gcb1[p->level][qq].cs, p->gcb1[p->level][qq].bc, gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,1);

    if(p->periodic1==1)
        gc_periodic(f, 1);

    if(p->periodic2==1)
        gc_periodic(f, 2);

    if(p->periodic3==1)
        gc_periodic(f, 3);

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol1(p,f);

    if(do_comms)
        gcparacox(p,f);
    #endif
}

void ghostcell::start2(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    #if USE_AMREX
    starttime=timer();
    f.FillDomainBoundary(gcv);
    endtime=timer();
    p->xtime+=endtime-starttime;
    p->gctime+=endtime-starttime;
    #else
    if(do_comms)
    {
        starttime=timer();
        gcparax(p,f,2);
        gcparacox(p,f);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    if(p->j_dir==1)
    {
        starttime=timer();
        QQGC2LOOP
            gcdistro2(f, p->gcb2[p->level][qq].i, p->gcb2[p->level][qq].j, p->gcb2[p->level][qq].k, p->gcb2[p->level][qq].cs, p->gcb2[p->level][qq].bc, gcv);
        endtime=timer();
        p->gctime+=endtime-starttime;

        // periodic ghostcells
        gcperiodicx(p,f,2);

        if(p->periodic1==1)
            gc_periodic(f, 1);

        if(p->periodic2==1)
            gc_periodic(f, 2);

        if(p->periodic3==1)
            gc_periodic(f, 3);
    }

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol2(p,f);

    if(do_comms)
        gcparacox(p,f);
    #endif
}

void ghostcell::start3(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    #if USE_AMREX
    starttime=timer();
    f.FillDomainBoundary(gcv);
    endtime=timer();
    p->xtime+=endtime-starttime;
    p->gctime+=endtime-starttime;
    #else
    if(do_comms)
    {
        starttime=timer();
        gcparax(p,f,3);
        gcparacox(p,f);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    starttime=timer();
    QQGC3LOOP
    if((p->gcb1[p->level][qq].cs!=2 && p->gcb1[p->level][qq].cs!=3) || p->j_dir==1)
        gcdistro3(f, p->gcb3[p->level][qq].i, p->gcb3[p->level][qq].j, p->gcb3[p->level][qq].k, p->gcb3[p->level][qq].cs, p->gcb3[p->level][qq].bc, gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,3);

    if(p->periodic1==1)
        gc_periodic(f, 1);

    if(p->periodic2==1)
        gc_periodic(f, 2);

    if(p->periodic3==1)
        gc_periodic(f, 3);

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol3(p,f);

    if(do_comms)
        gcparacox(p,f);
    #endif
}

void ghostcell::start4(lexer *p, field& f, int gcv, bool do_avgdown)
{
    //  MPI Boundary Swap
    #if USE_AMREX
    starttime=timer();
    f.FillDomainBoundary(gcv);
    endtime=timer();
    p->xtime+=endtime-starttime;
    // average_down overwrites covered coarse cells with the fine average. For the hydrostatic
    // pressure this mixes the fine roface/dz basis into the coarse column, breaking the coarse
    // grad(press) at the covered/non-covered surface boundary (the well-balancing seed). Callers
    // that need the coarse field kept self-consistent (press/press0) pass do_avgdown=false.
    if(do_avgdown)
    for(int lev=p->nlevs-2; lev>=0; --lev)
    {
        f.average_down_level(p, lev);
    }
    p->gctime+=endtime-starttime;
    #else
    if(do_comms)
    {
        starttime=timer();
        gcparax(p,f,4);
        gcparacox(p,f);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    starttime=timer();
    QQGC4LOOP
    {
        GCB4_TILE(qq);

        if((p->gcb1[p->level][qq].cs!=2 && p->gcb1[p->level][qq].cs!=3) || p->j_dir==1)
            gcdistro4(f, p->gcb4[p->level][qq].i, p->gcb4[p->level][qq].j, p->gcb4[p->level][qq].k, p->gcb4[p->level][qq].cs, p->gcb4[p->level][qq].bc, gcv);
    }
    GC_TILE_RESET;
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,4);

    if(p->periodic1==1)
        gc_periodic(f, 1);

    if(p->periodic2==1)
        gc_periodic(f, 2);

    if(p->periodic3==1)
        gc_periodic(f, 3);

    if(p->Y40==1 || p->Y40==3)
        dgcpol4(p,f);

    if(do_comms)
        gcparacox(p,f);
    #endif
}

#if USE_AMREX
void ghostcell::startBatch(lexer* p,
                           amrex::Vector<amrex::MultiFab>& shared_mf,
                           int scomp,
                           std::initializer_list<std::pair<field_amrex*, int>> fields_and_gcvs)
{
    starttime=timer();
    field_amrex::FillBoundaryBatch(p, shared_mf, scomp,
                                   static_cast<int>(fields_and_gcvs.size()));
    p->xtime+=timer()-starttime;
    starttime=timer();
    field_amrex::FillDomainBoundaryBatch(p, shared_mf, scomp, fields_and_gcvs, m_d_bcrec_batch);
    p->gctime+=timer()-starttime;
}
#endif

void ghostcell::start4_sum(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
    if(do_comms)
    {
        starttime=timer();
        gcparax4_sum(p,f,5);
        endtime=timer();
        p->xtime+=endtime-starttime;
    }

    starttime=timer();
    QQGC4LOOP
    {
        GCB4_TILE(qq);

        gcdistro4(f, p->gcb4[p->level][qq].i, p->gcb4[p->level][qq].j, p->gcb4[p->level][qq].k, p->gcb4[p->level][qq].cs, p->gcb4[p->level][qq].bc, gcv);
    }
    GC_TILE_RESET;
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,4);

    if(p->periodic1==1)
        gc_periodic(f, 1);

    if(p->periodic2==1)
        gc_periodic(f, 2);

    if(p->periodic3==1)
        gc_periodic(f, 3);

    if(do_comms)
        gcparacox(p,f);
}
