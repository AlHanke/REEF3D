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

#include"fnpf_vtp_fsf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<cstdio>
#include<cstring>

fnpf_vtp_fsf::fnpf_vtp_fsf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }
	
	printcount=0;
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_FNPF_VTP_FSF",0777);
    
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
	
	
}

void fnpf_vtp_fsf::start(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{	
    print2D(p,c,pgc);
}

void fnpf_vtp_fsf::print2D(lexer *p, fdm_fnpf *c, ghostcell* pgc)
{	
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    pgc->gcsl_start4(p,c->test2D,1);
    
    SLICELOOP4
    {
    if(c->breaking(i,j)>=1)
    c->breaking_print(i,j)=double(c->breaking(i,j));
        
    if(c->breaking(i,j)==0)
    c->breaking_print(i,j)=0.0;   
    }
    
    pgc->gcsl_start4(p,c->breaking_print,50);
    
	if(p->mpirank==0)
    pvtu(p,c,pgc);
    
	name_iter(p,c,pgc);
    
    // offsets
    n=0;
	offset[n]=0;
	++n;
	
	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
	++n;
    
    // Fifsf
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	
	// depth
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // breaking
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // coastline
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // test
    if(p->P23==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
    
    // Hs
    if(p->P110==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
	
    // Points
    offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
    ++n;
	// Cells
    offset[n]=offset[n-1] + 4*p->polygon_sum*3+4;
    ++n;
    offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
    
    //----------------------------------------------------------------------------

    std::stringstream result;

    beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);
    
    n=0;
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"Fifsf\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"depth\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"coastline\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P23==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    if(p->P110==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"Hs\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    result<<"</PointData>\n";

    points(result,offset,n);

    polys(result,offset,n);

    ending(result);

    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);
    
    //----------------------------------------------------------------------------
	
    //  Velocities
    iin=4*(p->pointnum2D)*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
    k = p->knoz-1;
	
	ffn=float(c->U[FIJKp1]);
    
    if(k==-1 && j==-1)
	ffn=float(c->U[FIJp1Kp1]);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);


	ffn=float(c->V[FIJKp1]);
    
    if(k==-1 && j==-1)
	ffn=float(c->V[FIJp1Kp1]);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);


	ffn=float(c->W[FIJKp1]);
    
    if(k==-1 && j==-1)
	ffn=float(c->W[FIJp1Kp1]);
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Fifsf
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
        
	ffn=float(p->sl_ipol4(c->Fifsf));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Eta
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->eta));
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
    
    //  Breaking
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
        
	ffn=float(p->sl_ipol4(c->breaking_print));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Coastline
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->coastline));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  test
    if(p->P23==1)
    {
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->test2D));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    //  Hs
    if(p->P110==1)
    {
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(c->Hs));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
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

	ffn=float(p->nhf_ipol4eta(p->wet,c->eta,c->bed)+p->wd);
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
