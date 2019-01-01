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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::paraini(lexer* p, fdm* a, ghostcell* pgc)
{
    p->xcoormax=pgc->globalmax(p->xcoormax);
    p->xcoormin=pgc->globalmin(p->xcoormin);

    p->ycoormax=pgc->globalmax(p->ycoormax);
    p->ycoormin=pgc->globalmin(p->ycoormin);

    p->zcoormax=pgc->globalmax(p->zcoormax);
    p->zcoormin=pgc->globalmin(p->zcoormin);
}
