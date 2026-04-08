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

#include"momentum_forcing.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"FSI.h"

void momentum_forcing::momentum_forcing_start(fdm* a, lexer* p, ghostcell *pgc, sixdof* p6dof, fsi* pfsi,
                                              field &u, field &v, field &w, field &fx, field &fy, field &fz, int iter, double alpha, bool final)
{
    double starttime=pgc->timer();

    // Forcing
    fx.setVal(0.0);
    fy.setVal(0.0);
    fz.setVal(0.0);

    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);

    pgc->solid_forcing(p,a,alpha,u,v,w,fx,fy,fz);

    p6dof->start_cfd(p,a,pgc,iter,u,v,w,fx,fy,fz,final);

    pfsi->forcing(p,a,pgc,alpha,u,v,w,fx,fy,fz,final);

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

    p->fbtime+=pgc->timer()-starttime;


    // ghostcell update
    pgc->solid_forcing_flag_update(p,a);
    pgc->gcdf_update(p,a);
    pgc->gcb_velflagio(p);
}
