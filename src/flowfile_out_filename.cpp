/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"flowfile_out.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void flowfile_out::filename(lexer *p, fdm *a, ghostcell *pgc)
{
	int num;
	
    num = p->count;

    // open file
    sprintf(name,"./REEF3D_FlowFile/REEF3D-flowfile-%0*i-%0*i.r3d",p->padding_time,n+1,p->padding_mpi,num);
  
}


