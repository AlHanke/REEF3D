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

#include"print_porous.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void print_porous::print_vtp(lexer *p, fdm *a, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_Porous",0777);
    
	sprintf(name,"./REEF3D_CFD_Porous/REEF3D_Porous-Object.vtp");

	ofstream result;
	result.open(name, ios::binary);

    n=0;

	offset[n]=0;
	++n;

    offset[n]=offset[n-1]+4*(vertice_num)*3 + 4;
    ++n;
    offset[n]=offset[n-1]+4*polygon_sum + 4;
    ++n;
    offset[n]=offset[n-1]+4*polygon_num + 4;
    ++n;
	//---------------------------------------------

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<PolyData>\n";
	result<<"<Piece NumberOfPoints=\""<<vertice_num<<"\" NumberOfPolys=\""<<polygon_num<<"\">\n";

    n=0;
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	++n;
    result<<"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";

	result<<"</Polys>\n";

    result<<"</Piece>\n";
    result<<"</PolyData>\n";

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";


//  XYZ
	iin=4*vertice_num*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<vertice_num;++n)
	{
	ffn=vertice[n][0];
	result.write((char*)&ffn, sizeof (float));

	ffn=vertice[n][1];
	result.write((char*)&ffn, sizeof (float));

	ffn=vertice[n][2];
	result.write((char*)&ffn, sizeof (float));
	}

//  Connectivity POLYGON
    iin=4*polygon_sum;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<polygon_num;++n)
	for(q=0;q<numvert[n];++q)
	{
	iin=polygon[n][q];
	result.write((char*)&iin, sizeof (int));
	}

//  Offset of Connectivity
    iin=4*polygon_num;
    result.write((char*)&iin, sizeof (int));
	iin=0;
	for(n=0;n<polygon_num;++n)
	{
	iin+=+ numvert[n];//polygon_offset[n];
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*polygon_num;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<polygon_num;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>\n";
    result<<"</VTKFile>\n";

	result.close();	
	
}
