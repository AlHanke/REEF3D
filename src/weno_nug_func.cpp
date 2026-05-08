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

        qfx.resize(i_size);
        qfy.resize(j_size);
        qfz.resize(k_size);

        cfx.resize(i_size);
        cfy.resize(j_size);
        cfz.resize(k_size);

        isfx.resize(i_size);
        isfy.resize(j_size);
        isfz.resize(k_size);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag = true;
    }
}

#if USE_AMREX
void weno_nug_func::rebuild_levels(lexer* p, int new_nlevs)
{
    if(!iniflag)
    {
        qfx.clear(); qfy.clear(); qfz.clear();
        cfx.clear(); cfy.clear(); cfz.clear();
        isfx.clear(); isfy.clear(); isfz.clear();

        i_size = max_i * new_nlevs;
        j_size = max_j * new_nlevs;
        k_size = max_k * new_nlevs;

        qfx.resize(i_size);
        qfy.resize(j_size);
        qfz.resize(k_size);

        cfx.resize(i_size);
        cfy.resize(j_size);
        cfz.resize(k_size);

        isfx.resize(i_size);
        isfy.resize(j_size);
        isfz.resize(k_size);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag=true;
    }
}
#endif
