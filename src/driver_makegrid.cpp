/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void driver::makegrid(lexer *p, ghostcell *pgc)
{	
    pgc->flagx(p,p->flagsf1);
    pgc->flagx(p,p->flagsf2);
    pgc->flagx(p,p->flagsf3);
    pgc->flagx(p,p->flagsf4);
    
    pgc->flagx(p,p->flag1);
    pgc->flagx(p,p->flag2);
    pgc->flagx(p,p->flag3);
    pgc->flagx(p,p->flag4);
    pgc->flagx(p,p->flag);
	pgc->gcxupdate(p);
    
    p->vecsize(pgc);
}
	
void driver::makegrid_cds()
{	
	pgc->sizeM_update(p,a);
}
	
