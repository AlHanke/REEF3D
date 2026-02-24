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

#include"weno_nug_func.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"

weno_nug_func::weno_nug_func(lexer* pp):epsilon(0.0),psi(1.0e-6)
{
    ini(pp);
    
    p=pp;
}

weno_nug_func::~weno_nug_func()
{
}

void weno_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        const int nlev_mult =
        #if USE_AMREX
            p->nlevs;
        #else
            1;
        #endif

        p->Darray(qfx,max_i*nlev_mult,2,6,2);
        p->Darray(qfy,max_j*nlev_mult,2,6,2);
        p->Darray(qfz,max_k*nlev_mult,2,6,2);

        p->Darray(cfx,max_i*nlev_mult,2,6);
        p->Darray(cfy,max_j*nlev_mult,2,6);
        p->Darray(cfz,max_k*nlev_mult,2,6);

        p->Darray(isfx,max_i*nlev_mult,2,6,3);
        p->Darray(isfy,max_j*nlev_mult,2,6,3);
        p->Darray(isfz,max_k*nlev_mult,2,6,3);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag=true;
    }
}

double ****weno_nug_func::qfx,****weno_nug_func::qfy,****weno_nug_func::qfz;
double ***weno_nug_func::cfx,***weno_nug_func::cfy,***weno_nug_func::cfz;
double ****weno_nug_func::isfx,****weno_nug_func::isfy,****weno_nug_func::isfz;
bool weno_nug_func::iniflag(false);
