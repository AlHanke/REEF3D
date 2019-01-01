/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"wallin_IM2.h"
#include"ikepsilon.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"convection.h"

wallin_IM2::wallin_IM2(lexer* p, fdm* a, ghostcell *pgc) : kepsilon_IM2(p,a,pgc),wallin(p,a),cmu(p->cmu)
{
	gcval_earsm=25;
}

wallin_IM2::~wallin_IM2()
{
}
void wallin_IM2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow)
{
	kepsilon_IM2::start(a,p,pconvec,pdiff,psolv,pgc,pflow);

	LOOP
	{
	tk = tau(a);
	sq(p,a);
	invar();
	beta(a);
	terms();
	aniso(a);
	}

	pgc->start4(p,rs11,gcval_earsm);
	pgc->start4(p,rs22,gcval_earsm);
	pgc->start4(p,rs33,gcval_earsm);
	pgc->start4(p,rs12,gcval_earsm);
	pgc->start4(p,rs13,gcval_earsm);
	pgc->start4(p,rs23,gcval_earsm);
}

void wallin_IM2::aniso(fdm* a)
{
	rs11(i,j,k) = (T1_11 + T3_11 + T4_11 + T6_11 + T9_11 )*kin(i,j,k);
	rs22(i,j,k) = (T1_22 + T3_22 + T4_22 + T6_22 + T9_22 )*kin(i,j,k);
	rs33(i,j,k) = -rs11(i,j,k)-rs22(i,j,k);
	rs12(i,j,k) = (T1_12 + T3_12 + T4_12 + T6_12 + T9_12 )*kin(i,j,k);
	rs13(i,j,k) = (T1_13 + T3_13 + T4_13 + T6_13 + T9_13 )*kin(i,j,k);
	rs23(i,j,k) = (T1_23 + T3_23 + T4_23 + T6_23 + T9_23 )*kin(i,j,k);
}

double wallin_IM2::tau(fdm* a)
{
		double t2;

       t = 6.0*sqrt(fabs(a->visc(i,j,k)/(fabs(eps(i,j,k))>(0.0)?(fabs(eps(i,j,k))):(1.0e20))));

		t2=(kin(i,j,k))/(fabs(eps(i,j,k))>(0.0)?(fabs(eps(i,j,k))):(1.0e20));

		if(t2>t)
		t=t2;

	return t;
}




