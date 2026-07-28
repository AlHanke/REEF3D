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

#include "geotopo.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "ioflow.h"
#include "reinitopo.h"

void geotopo::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow, reinitopo* preto, vrans* pvrans)
{
    if(p->toporead>0)
    {
        p->level = 0;
        // TILE_LOOP // IJK is not tile aware.
        IJKLOOP
        PBASECHECK
        a->topo(i,j,k) = p->flag_topo[IJK];
    }

    p->flag_topo.reset();

    if(p->S57>-1.0e20)
    {
        ALOOP
        a->topo(i,j,k)=-p->S57+p->ZP[KP];
    }

    preto->start(p,a,pgc,a->topo);

    if(p->S10==2)
    pflow->vrans_sed_update(p,a,pgc,pvrans);
}
