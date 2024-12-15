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

#include"fnpf_vtp_bed.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_vtp_bed::pvtp(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{	
	int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(name,"./REEF3D_FNPF_VTP_BED/REEF3D-FNPF-BED-%08i.pvtp",num);
	

	ofstream result;
	result.open(name);

	beginningParallel(p,result);
	
	result<<"<PPointData>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"depth\"/>\n";
	result<<"</PPointData>\n";
	
	pointsParallel(result);

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,c,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	endingParallel(result);

	result.close();

}

void fnpf_vtp_bed::piecename(lexer *p, fdm_fnpf *c, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-FNPF-BED-%08i-%06i.vtp",num,n+1);


}
