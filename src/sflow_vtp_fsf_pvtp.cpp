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
	
	sprintf(name,"./REEF3D_SFLOW_VTP_FSF/REEF3D-SFLOW-FSF-%08i.pvtp",num);

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<PPolyData  GhostLevel=\"0\">\n";
    
    if(p->P16==1)
    {
    result<<"<FieldData>\n";
    result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>\n";
    result<<"</FieldData>\n";
    }
	
	result<<"<PPoints>\n";
	result<<"<PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
	result<<"</PPoints>\n";
	
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
	
	result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\"/>\n";
	result<<"<DataArray type=\"Int32\" Name=\"offsets\" />\n";
    result<<"<DataArray type=\"Int32\" Name=\"types\" />\n";
	result<<"</Polys>\n";

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,b,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	result<<"</PPolyData>\n";
	result<<"</VTKFile>\n";

	result.close();

}

void sflow_vtp_fsf::piecename(lexer *p, fdm2D *b, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-SFLOW-FSF-%08i-%06i.vtp",num,n+1);
}
