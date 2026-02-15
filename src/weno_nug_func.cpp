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

weno_nug_func::weno_nug_func(lexer* p):epsilon(0.0),psi(1.0e-6)
{
    ini(p);

    weno_nug_func::p=p;
}

void weno_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        qfx.resize(max_i);
        qfy.resize(max_j);
        qfz.resize(max_k);

        cfx.resize(max_i);
        cfy.resize(max_j);
        cfz.resize(max_k);

        isfx.resize(max_i);
        isfy.resize(max_j);
        isfz.resize(max_k);

        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag = true;
    }
}
