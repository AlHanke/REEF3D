
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"sflow_state.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void sflow_state::ini_mainheader(lexer *p, fdm2D *b, ghostcell *pgc)
{
    ofstream mainout;
    
    // open file
	if(p->P14==0)
    mainout.open("REEF3D-SFLOW_State_Mainheader.r3d", ios::binary);

	if(p->P14==1)
	mainout.open("./REEF3D_SFLOW_STATE/REEF3D-SFLOW_State_Mainheader.r3d", ios::binary);


    // ini write

    // flag: is process within P43 bounds
    for(int qn=0;qn<p->M10;++qn)
    {
    iin = flag_all[qn];
    mainout.write((char*)&iin, sizeof (int));
    }

    mainout.close();
}

void sflow_state::write_mainheader(lexer *p, fdm2D *c, ghostcell *pgc)
{
    ofstream mainout;
    
    // open file
	if(p->P14==0)
    mainout.open("REEF3D-SFLOW_State_Mainheader.r3d", ios::binary | ios::app);

	if(p->P14==1)
	mainout.open("./REEF3D_SFLOW_STATE/REEF3D-SFLOW_State_Mainheader.r3d", ios::binary | ios::app);

    iin=p->count;
    mainout.write((char*)&iin, sizeof (int));

	ddn=p->simtime;
    mainout.write((char*)&ddn, sizeof (double));

    mainout.close();
}