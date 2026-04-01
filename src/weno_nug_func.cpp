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

        qfx.resize(max_i*nlev_mult);
        qfy.resize(max_j*nlev_mult);
        qfz.resize(max_k*nlev_mult);

        cfx.resize(max_i*nlev_mult);
        cfy.resize(max_j*nlev_mult);
        cfz.resize(max_k*nlev_mult);

        isfx.resize(max_i*nlev_mult);
        isfy.resize(max_j*nlev_mult);
        isfz.resize(max_k*nlev_mult);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag = true;
    }
}
