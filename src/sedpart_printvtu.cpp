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

#include"sedpart.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sedpart::print_particles(lexer* p, fdm* a, ghostcell* pgc)
{
    if(((p->count%p->Q181==0 && p->Q182<0.0 && p->Q180==1 )|| (p->count==0 &&  p->Q182<0.0 && p->Q180==1)) && p->Q181>0)
	{
    print_vtu(p,a,pgc);
	++printcount;
	}
    
    if((p->simtime>p->fsfprinttime && p->Q182>0.0 && p->Q180==1) || (p->count==0 &&  p->Q182>0.0))
    {
    print_vtu(p,a,pgc);
    p->partprinttime+=p->Q182;
    }
    
}

void sedpart::print_vtu(lexer* p, fdm* a, ghostcell* pgc)
{
	int numpt=0;
	const int print_flag=p->Q183;

	if(print_flag==0)
		numpt=PP.size;
	else
		PARTLOOP
			if(PP.Flag[n]>print_flag)
				numpt++;

	cout<<"PSedACTIVE-"<<p->mpirank<<": "<<numpt<<"|"<<PP.capacity<<endl;

	int count;
	int n=0;
    int offset[100];
	int iin;
	float ffn;
	
	if(p->mpirank==0)
	pvtu_pos(p,a,pgc);


    header_pos(p,a,pgc);

	ofstream result;
	result.open(name, ios::binary);

    

	offset[n]=0;
	++n;
	
	//offset[n]=offset[n-1]+4*(numpt)+4; //lsv
	//++n;
	// offset[n]=offset[n-1]+4*(numpt)+4; //radius
	// ++n;
	//offset[n]=offset[n-1]+4*(numpt)+4; //correction
	//++n;	
	
	// end scalars
    offset[n]=offset[n-1]+4*(numpt)*3+4; //xyz
    ++n;
    offset[n]=offset[n-1]+4*(numpt)*2+4; //connectivitey
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //offset connectivity
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //cell type
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfCells=\""<<numpt<<"\">"<<endl;
	
	
	// result<<"<PointData >"<<endl;
	// //result<<"<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    // //++n;
    // // result<<"<DataArray type=\"Float32\" Name=\"radius\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    // // ++n;
	// //result<<"<DataArray type=\"Float32\" Name=\"correction\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    // //++n;
	// result<<"</PointData>"<<endl;
	
	

    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
	

    result<<"<Cells>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</UnstructuredGrid>"<<endl;

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	

//  lsv
    // iin=4*(numpt);
    // result.write((char*)&iin, sizeof (int));
	// PARTLOOP
    // if(PP.Flag[n]>0)
	// {
	// ffn=float(f[n][3]);
	// result.write((char*)&ffn, sizeof (float));
	// }
	
//  radius
    // iin=4*(numpt);
    // result.write((char*)&iin, sizeof (int));
	// PARTLOOP
	// 	if(PP.Flag[n]>0)
	// 	{
	// 		ffn=float(1);
	// 		result.write((char*)&ffn, sizeof (float));
	// 	}

//  correction
    // iin=4*(numpt);
    // result.write((char*)&iin, sizeof (int));
	// PARTLOOP
    // if(PP.Flag[n]>0)
	// {
	// 	if(sign==1)
	// 	{
	// 	if(f[n][3]<=-f[n][4])
	// 	ffn=float(1.0);
		
	// 	if(f[n][3]>-f[n][4])
	// 	ffn=float(0.0);
	// 	}
		
	// 	if(sign==2)
	// 	{
	// 	if(f[n][3]>=f[n][4])
	// 	ffn=float(1.0);
		
	// 	if(f[n][3]<f[n][4])
	// 	ffn=float(0.0);
	// 	}
		
	// result.write((char*)&ffn, sizeof (float));
	// }

//  XYZ
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    PARTLOOP
    if(PP.Flag[n]>print_flag)
	{
	ffn=float(PP.X[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.Y[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.Z[n]);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Connectivity
	count=0;
    iin=4*(numpt)*2;
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
	if(PP.Flag[n]>print_flag)
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
	PARTLOOP
    if(PP.Flag[n]>print_flag)
	{
	iin=(count+1)*2;
	result.write((char*)&iin, sizeof (int));
	++count;
	}


//  Cell types
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
    if(PP.Flag[n]>print_flag)
	{
	iin=1;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
	}

