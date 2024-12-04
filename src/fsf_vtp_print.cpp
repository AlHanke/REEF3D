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

#include"fsf_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sstream>
#include<cstdio>
#include<cstring>

void fsf_vtp::print(lexer* p, fdm* a, ghostcell *pgc)
{
	if(p->mpirank==0)
    pvtp(p,a,pgc);
	
	name_iter(p,a,pgc);

	//---------------------------------------------
	
	polygon_num=facount;
	
	polygon_sum=0;
	for(n=0;n<polygon_num;++n)
	polygon_sum+=numpt[n];
	
	point_num = ccptcount;
	
	//---------------------------------------------
    n=0;
	offset[n]=0;
	++n;
	//Data
	offset[n]=offset[n-1] + 4*point_num*3+ 4;
    ++n;
	offset[n]=offset[n-1] + 4*point_num+ 4;
    ++n;
	//End Data
    // Points
    offset[n]=offset[n-1] + 4*point_num*3 + 4;
    ++n;
    // Polys connectivity
    offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
    // Polys offset
    offset[n]=offset[n-1] + 4*polygon_num+ 4;
    ++n;
	//---------------------------------------------
	
	std::stringstream result;

    beginning(p,result,point_num,0,0,0,polygon_num);

    n=0;
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    points(result,offset,n);

    polys(result,offset,n);
	
    ending(result);

    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);

    //----------------------------------------------------------------------------
	
    //  Velocity
	iin=4*point_num*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<point_num;++n)
	{
	ffn=float(p->ccipol1(a->u,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->ccipol2(a->v,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->ccipol3(a->w,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
	
	
    //  Elevation
	iin=4*point_num;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<point_num;++n)
	{
	ffn=ccpt[n][2];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  XYZ
	iin=4*point_num*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<point_num;++n)
	{
	ffn=ccpt[n][0];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=ccpt[n][1];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=ccpt[n][2];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  Connectivity POLYGON
    iin=4*polygon_sum;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<polygon_num;++n)
	{
		if(numpt[n]==3)
		{
		iin=facet[n][0];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		
		iin=facet[n][1];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		
		iin=facet[n][2];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		}
		
		if(numpt[n]==4)
		{
		iin=facet[n][0];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		
		iin=facet[n][1];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		
		iin=facet[n][3];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		
		iin=facet[n][2];
		std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
		}
	}

    //  Offset of Connectivity
    iin=4*polygon_num;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	iin=0;
	for(n=0;n<polygon_num;++n)
	{
	iin+= numpt[n];
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	}

	footer(buffer,m);

	// Open File
    FILE* file = std::fopen(name, "w");
    std::fwrite(buffer.data(), buffer.size(), 1, file);
    std::fclose(file);
}
