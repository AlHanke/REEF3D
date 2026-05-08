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

weno_nug_func::weno_nug_func(lexer* p)
{
    ini(p);

    weno_nug_func::p=p;
    #if USE_AMREX
    p->register_weno5(this);
    #endif
}

weno_nug_func::~weno_nug_func()
{
    #if USE_AMREX
    p->deregister_weno5(this);
    #endif
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

        i_size = max_i * nlev_mult;
        j_size = max_j * nlev_mult;
        k_size = max_k * nlev_mult;

        p->Darray(qfx,i_size,2,6,2);
        p->Darray(qfy,j_size,2,6,2);
        p->Darray(qfz,k_size,2,6,2);

        p->Darray(cfx,i_size,2,6);
        p->Darray(cfy,j_size,2,6);
        p->Darray(cfz,k_size,2,6);

        p->Darray(isfx,i_size,2,6,3);
        p->Darray(isfy,j_size,2,6,3);
        p->Darray(isfz,k_size,2,6,3);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag=true;
    }
}

#if USE_AMREX
void weno_nug_func::rebuild_levels(lexer* p, int new_nlevs)
{
    if(!iniflag)
    {
        p->del_Darray(qfx,i_size,2,6,2); p->del_Darray(qfy,j_size,2,6,2); p->del_Darray(qfz,k_size,2,6,2);
        p->del_Darray(cfx,i_size,2,6); p->del_Darray(cfy,j_size,2,6); p->del_Darray(cfz,k_size,2,6);
        p->del_Darray(isfx,i_size,2,6,3); p->del_Darray(isfy,j_size,2,6,3); p->del_Darray(isfz,k_size,2,6,3);

        i_size = max_i * new_nlevs;
        j_size = max_j * new_nlevs;
        k_size = max_k * new_nlevs;

        p->Darray(qfx, i_size, 2, 6, 2);
        p->Darray(qfy, j_size, 2, 6, 2);
        p->Darray(qfz, k_size, 2, 6, 2);

        p->Darray(cfx, i_size, 2, 6);
        p->Darray(cfy, j_size, 2, 6);
        p->Darray(cfz, k_size, 2, 6);

        p->Darray(isfx, i_size, 2, 6, 3);
        p->Darray(isfy, j_size, 2, 6, 3);
        p->Darray(isfz, k_size, 2, 6, 3);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag=true;
    }
}
#endif

double ****weno_nug_func::qfx,****weno_nug_func::qfy,****weno_nug_func::qfz;
double ***weno_nug_func::cfx,***weno_nug_func::cfy,***weno_nug_func::cfz;
double ****weno_nug_func::isfx,****weno_nug_func::isfy,****weno_nug_func::isfz;
bool weno_nug_func::iniflag(false);
