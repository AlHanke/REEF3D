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

#include"nhflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<cstdio>
#include<cstring>

nhflow_vtp_fsf::nhflow_vtp_fsf(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }
	
	printcount=0;
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_VTP_FSF",0777);
    
    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
    
    if(p->P131==1)
    {
    p->Iarray(wetmax,p->imax*p->jmax);
    
    SLICELOOP4
    wetmax[IJ] = 0;
    
    pgc->gcsl_start4Vint(p,wetmax,50);
    }
}

nhflow_vtp_fsf::~nhflow_vtp_fsf()
{
}

void nhflow_vtp_fsf::start(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{	
    print2D(p,d,pgc,psed);
}

void nhflow_vtp_fsf::print2D(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{	
    //pgc->gcsl_start4(p,d->eta,gcval_eta);
    //pgc->gcsl_start4(p,d->Fifsf,gcval_fifsf);
    
    pgc->gcsl_start4(p,d->test2D,1);
    
    SLICELOOP4
    {
    if(d->breaking(i,j)>=1)
    d->breaking_print(i,j)=double(d->breaking(i,j));
        
    if(d->breaking(i,j)==0)
    d->breaking_print(i,j)=0.0;   
    }
    
    //pgd->gcsl_start4(p,d->breaking_print,50);
    

	if(p->mpirank==0)
    pvtu(p,d,pgc,psed);

	name_iter(p,d,pgc);
    
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
	
	// WL
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // breaking
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // coastline
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // wetdry
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    
    // test
    if(p->P23==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }
    
    // fb
    if(p->P28==1)
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
    
    // wetdry_max
    if(p->P131==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }

    // Points
    offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
    ++n;
	
	// Polys
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
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"WL\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"coastline\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"wetdry\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    
    if(p->P23==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    
    if(p->P28==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"fb\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    
    if(p->P110==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"Hs\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    }
    
    if(p->P131==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"wetdry_max\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
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
    
	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->U[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->U[IJK]+d->U[IJp1K]));
    
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);


	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->V[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->V[IJK]+d->V[IJp1K]));
    
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);


	if(p->j_dir==0)
    {
	jj=j;
    j=0;
	ffn=float(d->W[IJK]);
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.5*(d->W[IJK]+d->W[IJp1K]));
    
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Eta
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4eta_wd(p->wet,d->eta,d->bed));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  WL
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->WL));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Breaking
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
        
	ffn=float(p->sl_ipol4(d->breaking_print));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Coastline
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->coastline));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Wetdry
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
    ffn = 0.25*float((p->wet[IJ]+p->wet[Ip1J]+p->wet[IJp1]+p->wet[Ip1Jp1]));
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
	ffn=float(p->sl_ipol4(d->test2D));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    //  fb
    if(p->P28==1)
    {
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(d->fs));
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
	ffn=float(p->sl_ipol4(d->Hs));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    //  wetdry_max
    if(p->P131==1)
    {
	iin=4*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
    ffn = 0.25*float((wetmax[IJ]+wetmax[Ip1J]+wetmax[IJp1]+wetmax[Ip1Jp1]));
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
	ffn=p->XN[IP1];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=p->YN[JP1];
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

    //ddn=float(p->sl_ipol4(d->eta) + p->wd);
    
    //ddn=float(0.5*(d->eta(i,j) + d->eta(i,j))  + p->wd);
    ffn=p->nhf_ipol4eta(p->wet,d->eta, d->bed)+p->wd;
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
	iin=int(d->nodeval2D(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(d->nodeval2D(i,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(d->nodeval2D(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	
	// Triangle 2
	iin=int(d->nodeval2D(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(d->nodeval2D(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(d->nodeval2D(i-1,j))-1;
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

void nhflow_vtp_fsf::preproc(lexer *p, fdm_nhf *d, ghostcell* pgc)
{	
    if(p->P131==1)
    {
    SLICELOOP4
    wetmax[IJ] = MAX(wetmax[IJ],p->wet[IJ]);
    
    pgc->gcsl_start4Vint(p,wetmax,50);    
    }
    
}


