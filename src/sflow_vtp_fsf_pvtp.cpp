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

#include"sflow_vtp_fsf.h"
#include"sflow_turbulence.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment.h"

void sflow_vtp_fsf::pvtp(lexer *p, fdm2D* b, ghostcell* pgc, sflow_turbulence *pturb, sediment *psed)
{	
	int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;
	
	snprintf(name,sizeof(name),"./REEF3D_SFLOW_VTP_FSF/REEF3D-SFLOW-FSF-%08i.pvtp",num);

	ofstream result;
	result.open(name);

	beginningParallel(p,result);
	
	result<<"<PPointData>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>\n";
    pturb->name_pvtp(p,b,pgc,result);
	result<<"<PDataArray type=\"Float32\" Name=\"eta\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"waterlevel\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"wetdry\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"breaking\"/>\n";
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";
    if(p->P28==1)
    result<<"<PDataArray type=\"Float32\" Name=\"fb\"/>\n";
    if(p->P110==1)
    result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>\n";
	result<<"</PPointData>\n";
	
	pointsParallel(result);

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,b,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	endingParallel(result);

	result.close();

}

void sflow_vtp_fsf::piecename(lexer *p, fdm2D *b, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	snprintf(pname,sizeof(pname),"REEF3D-SFLOW-FSF-%08i-%06i.vtp",num,n+1);
}
