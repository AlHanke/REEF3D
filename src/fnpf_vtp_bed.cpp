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
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<cstdio>
#include<cstring>

fnpf_vtp_bed::fnpf_vtp_bed(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }
	
	printcount=0;
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_VTP_BED",0777);
}

void fnpf_vtp_bed::start(lexer *p, fdm_fnpf *c, ghostcell* pgc, ioflow *pflow)
{	
    print2D(p,c,pgc);
}

void fnpf_vtp_bed::print2D(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{	
    
	if(p->mpirank==0)
    pvtp(p,c,pgc);
    
	name_iter(p,c,pgc);
    
    // offsets
    n=0;
	offset[n]=0;
	++n;
    
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	
	// depth
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // Points
    offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
    ++n;
	// Polys
    offset[n]=offset[n-1] + 4*p->polygon_sum*3+4;
    ++n;
    offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;

    // Open File
	std::stringstream result;
    
    beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);
    
    n=0;
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    points(result,offset,n);

    polys(result,offset,n);

    ending(result);
    
    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);
    
    //----------------------------------------------------------------------------
	
    //  Elevation
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->bed));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
	
	//  Depth
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->depth));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  XYZ
	iin=4*(p->pointnum2D)*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
    
	ffn=float(p->XN[IP1]);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->YN[JP1]);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->sl_ipol4(c->bed));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  Connectivity
    iin=4*(p->polygon_sum)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(c->nodeval2D(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(c->nodeval2D(i,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(c->nodeval2D(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	
	// Triangle 2
	iin=int(c->nodeval2D(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(c->nodeval2D(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(c->nodeval2D(i-1,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	}
    
    
    //  Offset of Connectivity
    iin=4*(p->polygon_sum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=(n+1)*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	}

    footer(buffer,m);

	// Open File
    FILE* file = std::fopen(name, "w");
    std::fwrite(buffer.data(), buffer.size(), 1, file);
    std::fclose(file);
	
	++printcount;

}
