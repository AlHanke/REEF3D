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

#include "momentum_forcing.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "6DOF.h"
#include "FSI.h"
#include "field1.h"
#include "field2.h"
#include "field3.h"
#include <algorithm>

void momentum_forcing::momentum_forcing_start(lexer* p, fdm* a, ghostcell *pgc, sixdof* p6dof, fsi* pfsi,
                                              field1 &u, field2 &v, field3 &w, field1 &fx, field2 &fy, field3 &fz, int iter, double alpha, bool final)
{
    double starttime=pgc->timer();

    // Forcing
    fx.setVal(0.0, true);
    fy.setVal(0.0, true);
    fz.setVal(0.0, true);

    double setzero_time = pgc->timer();

    pgc->solid_forcing(p,a,alpha,u,v,w,fx,fy,fz);

    double gc_forcing_time = pgc->timer();

    p6dof->start_cfd(p,a,pgc,iter,u,v,w,fx,fy,fz,final);

    double sixdof_time = pgc->timer();

    pfsi->forcing(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);

    double compute_time = pgc->timer();

    ULOOP
    {
        const double val = alpha*CPOR1*fx(i,j,k);

        u(i,j,k) += val*p->dt;

        if(p->count<10)
        a->maxF = MAX(fabs(val), a->maxF);

        p->fbmax = MAX(fabs(val), p->fbmax);
    }

    VLOOP
    {
        const double val = alpha*CPOR2*fy(i,j,k);

        v(i,j,k) += val*p->dt;

        if(p->count<10)
        a->maxG = MAX(fabs(val), a->maxG);

        p->fbmax = MAX(fabs(val), p->fbmax);
    }

    WLOOP
    {
        const double val = alpha*CPOR3*fz(i,j,k);

        w(i,j,k) += val*p->dt;

        if(p->count<10)
        a->maxH = MAX(fabs(val), a->maxH);

        p->fbmax = MAX(fabs(val), p->fbmax);
    }

    double apply_time = pgc->timer();

    p->fbtime+=pgc->timer()-starttime;

    double gc_update_time2 = apply_time;

    // ghostcell update
    if(iter==0)
    {
        pgc->solid_forcing_flag_update(p,a);

        gc_update_time2 = pgc->timer();

        if(p->S10>0 || p->X10>0 || p->Z10>0)
        pgc->gcdf_update(p,a);
    }

    double gc_update_time = pgc->timer();

    if(p->mpirank==0 && p->count%p->P12==0 && p->count>0 && std::getenv("REEF_timing"))
    {
        const int precision = 5;
        const double total = gc_update_time - starttime;
        const double zero = setzero_time-starttime;
        const double forcing = gc_forcing_time-setzero_time;
        const double compute = compute_time-sixdof_time;
        std::cout<<"\nTiming for momentum forcing iteration "<<iter<<"\n"
        <<"total time:     "<<std::setprecision(precision)<<total<<"\n"
        <<"\tset zero:     "<<std::setprecision(precision)<<zero<<":"<<100.0*(zero/total)<<"\n"
        <<"\tforcing:      "<<std::setprecision(precision)<<forcing<<":"<<100.0*forcing/total<<"\n"
        <<"\tsixdof:       "<<std::setprecision(precision)<<sixdof_time-gc_forcing_time<<":"<<100.0*(sixdof_time-gc_forcing_time)/total<<"\n"
        <<"\tcompute time: "<<std::setprecision(precision)<<compute<<":"<<100.0*compute/total<<"\n"
        <<"\tapply time:   "<<std::setprecision(precision)<<apply_time-compute_time<<":"<<100.0*(apply_time-compute_time)/total<<"\n"
        <<"\tgc_sf:        "<<std::setprecision(precision)<<gc_update_time2-apply_time<<":"<<100.0*(gc_update_time2-apply_time)/total<<"\n"
        <<"\tgc update:    "<<std::setprecision(precision)<<gc_update_time-gc_update_time2<<":"<<100.0*(gc_update_time-gc_update_time2)/total<<"\n"
        <<std::endl;
    }
}
