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

#include"weno3_nug_func.h"
#include"lexer.h"

weno3_nug_func::weno3_nug_func(lexer* p)
{
    ini(p);
    weno3_nug_func::p=p;
}

void weno3_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        const int nlev_mult =
        #if USE_AMREX
                p->nlevs;
        #else
                1;
        #endif

        p->Darray(qfx,max_i*nlev_mult,2,4,2);
        p->Darray(qfy,max_j*nlev_mult,2,4,2);
        p->Darray(qfz,max_k*nlev_mult,2,4,2);

        p->Darray(cfx,max_i*nlev_mult,2,4);
        p->Darray(cfy,max_j*nlev_mult,2,4);
        p->Darray(cfz,max_k*nlev_mult,2,4);

        p->Darray(isfx,max_i*nlev_mult,2,4);
        p->Darray(isfy,max_j*nlev_mult,2,4);
        p->Darray(isfz,max_k*nlev_mult,2,4);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag=true;
    }
}

double ****weno3_nug_func::qfx,****weno3_nug_func::qfy,****weno3_nug_func::qfz;
double ***weno3_nug_func::cfx,***weno3_nug_func::cfy,***weno3_nug_func::cfz;
double ***weno3_nug_func::isfx,***weno3_nug_func::isfy,***weno3_nug_func::isfz;
bool weno3_nug_func::iniflag(false);
