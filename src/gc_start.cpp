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

void ghostcell::start1(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
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
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
        gcdistro1(f, p->gcb1[qq][0], p->gcb1[qq][1], p->gcb1[qq][2], p->gcb1[qq][3], p->gcb1[qq][4], p->gcd1[qq], gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,1);

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol1(p,f);

    if(do_comms)
        gcparacox(p,f);
}

void ghostcell::start2(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
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
            gcdistro2(f, p->gcb2[qq][0], p->gcb2[qq][1], p->gcb2[qq][2], p->gcb2[qq][3], p->gcb2[qq][4], p->gcd2[qq], gcv);
        endtime=timer();
        p->gctime+=endtime-starttime;

        // periodic ghostcells
        gcperiodicx(p,f,2);
    }

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol2(p,f);

    if(do_comms)
        gcparacox(p,f);
}

void ghostcell::start3(lexer *p, field& f, int gcv)
{
    //  MPI Boundary Swap
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
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
        gcdistro3(f, p->gcb3[qq][0], p->gcb3[qq][1], p->gcb3[qq][2], p->gcb3[qq][3], p->gcb3[qq][4], p->gcd3[qq], gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,3);

    if((p->Y40==1 || p->Y40==3) && p->j_dir==1)
        dgcpol3(p,f);

    if(do_comms)
        gcparacox(p,f);
}

void ghostcell::start4(lexer *p, field &f, int gcv)
{
    //  MPI Boundary Swap
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
    if((p->gcb1[qq][3]!=2 && p->gcb1[qq][3]!=3) || p->j_dir==1)
        gcdistro4(f, p->gcb4[qq][0], p->gcb4[qq][1], p->gcb4[qq][2], p->gcb4[qq][3], p->gcb4[qq][4], p->gcd4[qq], gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,4);

    if(p->Y40==1 || p->Y40==3)
        dgcpol4(p,f);

    if(do_comms)
        gcparacox(p,f);
}

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
        gcdistro4(f, p->gcb4[qq][0], p->gcb4[qq][1], p->gcb4[qq][2], p->gcb4[qq][3], p->gcb4[qq][4], p->gcd4[qq], gcv);
    endtime=timer();
    p->gctime+=endtime-starttime;

    // periodic ghostcells
    gcperiodicx(p,f,4);

    if(do_comms)
        gcparacox(p,f);
}
