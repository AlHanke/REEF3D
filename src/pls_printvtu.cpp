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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_pls::print_vtu(lexer* p, fdm* a, ghostcell* pgc,double** f,int *flag,int active, int sign)
{
	int numpt=0;
	int count;
	
	for(n=0;n<active;++n)
    if(flag[n]>0)
	++numpt;
	
	
	
    if(p->count>0)
    
	
	if(p->mpirank==0)
	{
		if(sign==1)
		pvtu_pos(p,a,pgc);
		
		if(sign==2)
		pvtu_neg(p,a,pgc);	
	}
	
	if(sign==1)
    header_pos(p,a,pgc);
	
	if(sign==2)
    header_neg(p,a,pgc);


	ofstream result;
	result.open(name, ios::binary);

    n=0;

	offset[n]=0;
	++n;
	
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;	
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;	
	
	// end scalars
    offset[n]=offset[n-1]+4*(numpt)*3+4;
    ++n;
    offset[n]=offset[n-1]+4*(numpt)*2+4;
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<UnstructuredGrid>\n";
	result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfCells=\""<<numpt<<"\">\n";
	
	
	result<<"<PointData>\n";
	result<<"<DataArray type=\"Float32\" Name=\"phi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"correction\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"</PointData>\n";
	
	
	

    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";
	
	

    result<<"<Cells>\n";
	result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	++n;
    result<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"</Cells>\n";

    result<<"</Piece>\n";
    result<<"</UnstructuredGrid>\n";

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	

//  lsv
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][3]);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  radius
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][4]);
	result.write((char*)&ffn, sizeof (float));
	}

//  correction
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
		if(sign==1)
		{
		if(f[n][3]<=-f[n][4])
		ffn=float(1.0);
		
		if(f[n][3]>-f[n][4])
		ffn=float(0.0);
		}
		
		if(sign==2)
		{
		if(f[n][3]>=f[n][4])
		ffn=float(1.0);
		
		if(f[n][3]<f[n][4])
		ffn=float(0.0);
		}
		
	result.write((char*)&ffn, sizeof (float));
	}

//  XYZ
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][0]+p->originx);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(f[n][1]+p->originy);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(f[n][2]+p->originz);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Connectivity
	count=0;
    iin=4*(numpt)*2;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
	if(flag[n]>0)
	{
	iin=int(0);
	result.write((char*)&iin, sizeof (int));

	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

//  Offset of Connectivity
	count=0;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	iin=(count+1)*2;
	result.write((char*)&iin, sizeof (int));
	++count;
	}


//  Cell types
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	iin=1;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>\n";
    result<<"</VTKFile>\n";

	result.close();
}

