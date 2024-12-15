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
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment.h"
#include"sflow_print_wsf.h"
#include"sflow_print_wsf_theory.h"
#include"sflow_print_wsfline.h"
#include"sflow_print_wsfline_y.h"
#include"sflow_print_probe_da.h"
#include"sflow_print_bed.h"
#include"sflow_print_bedline.h"
#include"sflow_print_bedline_y.h"
#include"sflow_turbulence.h"
#include"sflow_state.h"
#include"fnpf_print_Hs.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<cstdio>
#include<cstring>
#include<chrono>

sflow_vtp_fsf::sflow_vtp_fsf(lexer *p, fdm2D *b, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }

	p->printcount=0;

	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_SFLOW_VTP_FSF",0777);


	pwsf=new sflow_print_wsf(p,b);

	pwsf_theory=new sflow_print_wsf_theory(p,b,pgc);

    pwsfline=new sflow_print_wsfline(p,b,pgc);

    pwsfline_y=new sflow_print_wsfline_y(p,b,pgc);

    pprobe=new sflow_print_probe_da(p,b,pgc);
    
    pbed=new sflow_print_bed(p,b);

    pbedline=new sflow_print_bedline(p,b,pgc);

    pbedline_y=new sflow_print_bedline_y(p,b,pgc);
    
    if(p->P40>0)
	pstate=new sflow_state(p,b,pgc);
    
    if(p->P110==1)
    phs = new fnpf_print_Hs(p,b->Hs);
}

sflow_vtp_fsf::~sflow_vtp_fsf()
{
}

void sflow_vtp_fsf::start(lexer *p, fdm2D* b, ghostcell* pgc, ioflow *pflow, sflow_turbulence *pturb, sediment *psed)
{
	// Print out based on iteration
    if((p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)  || (p->count==0 &&  p->P30<0.0))
    {
    print2D(p,b,pgc,pturb,psed);
    }

    // Print out based on time
    if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
    {
    print2D(p,b,pgc,pturb,psed);

    p->printtime+=p->P30;
    }

	// WSF Gages
    if(p->P51>0)
    pwsf->height_gauge(p,b,pgc,b->eta);

    if(p->P50>0)
    pwsf_theory->height_gauge(p,b,pgc,pflow);

    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,b,pgc,pflow,b->eta);

    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,b,pgc,pflow,b->eta);
    
    // DA Gages
    if(p->P63>0 && p->count%p->P54==0)
    pprobe->start(p,b,pgc);
    
    // Hs
    if(p->P110==1)
    phs->start(p,pgc,b->eta,b->Hs);
    
    // BED Gages
    if(((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47) ) && p->S10>0)
    if((p->S42==1 && p->count%p->S44==0 && p->sediter%p->P120==0) || (p->S42==2 && p->simtime>=p->sedsimtime && p->sediter%p->P120==0) || (p->S42==3  && p->simtime/p->wT>=p->sedwavetime && p->sediter%p->P120==0))
    {      
    if(p->P121>0)
    pbed->height_gauge(p,b,pgc,b->bed);

    if(p->P123>0)
    pbedline->start(p,b,pgc,pflow,b->bed);

    if(p->P124>0)
    pbedline_y->start(p,b,pgc,pflow,b->bed);
    }
    
    // Print state out based on iteration
    if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && (p->P46==0 || (p->count>=p->P46_is && p->count<<p->P46_ie)))
    {
    pstate->write(p,b,pgc);
    }

    // Print state out based on time
    if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0 && (p->P47==0 || (p->count>=p->P47_ts && p->count<<p->P47_te)))
    {
    pstate->write(p,b,pgc);

    p->stateprinttime+=p->P42;
    }
    
    if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
    p->probeprinttime+=p->P55;
}

void sflow_vtp_fsf::print2D(lexer *p, fdm2D* b, ghostcell* pgc, sflow_turbulence *pturb, sediment *psed)
{
	if(p->mpirank==0)
    pvtp(p,b,pgc,pturb,psed);

	name_iter(p,b,pgc);
    
    if(p->origin_i==0 && p->origin_j==0)
    {
    i=0;
    j=0;
	p->wet[Im1Jm1] = p->wet[IJ];
    }
    
    if(p->origin_i==0 && p->gknoy==p->knoy+p->origin_j)
    {
    i=0;
    j=p->knoy-1;
	p->wet[Im1Jp1] = p->wet[IJ];
    }
    
    if(p->gknox==p->knox+p->origin_i && p->origin_j==0)
    {
    i=p->knox-1;
    j=0;
	p->wet[Ip1Jm1] = p->wet[IJ];
    }
    
    if(p->gknox==p->knox+p->origin_i && p->gknoy==p->knoy+p->origin_j)
    {
    i=p->knox-1;
    j=p->knoy-1;
	p->wet[Ip1Jp1] = p->wet[IJ];
    }

    auto start = chrono::high_resolution_clock::now();

    // offsets
    n=0;
	offset[n]=0;
	++n;

	// velocity
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+4;
	++n;

    // pressure
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    
    // eddyv
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    
    // k and eps
	pturb->offset_ParaView_2D(p,offset,n);
    
    // eta
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    
    // waterlevel
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    
    // wetdry
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;

    // breaking
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;

    // test
    if(p->P23==1)
    {
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    }
    
    // fb
    if(p->P28==1)
    {
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    }
    
    // Hs
    if(p->P110==1)
	{
	offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)+4;
	++n;
    }

    // Points
    offset[n]=offset[n-1]+sizeof(float)*(p->pointnum2D)*3+4;
    ++n;
	// Polys
    // Connectivity
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum*3+4;
    ++n;
    // Offsets
    offset[n]=offset[n-1] + sizeof(int)*p->polygon_sum+4;
    ++n;

    //----------------------------------------------------------------------------
    std::stringstream result;
    beginning(p,result,p->pointnum2D,0,0,0,p->polygon_sum);

    n=0;
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    pturb->name_vtp(p,b,pgc,result,offset,n);
    result<<"<DataArray type=\"Float32\" Name=\"eta\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"waterlevel\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"wetdry\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
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
    result<<"</PointData>\n";

    points(result,offset,n);

    polys(result,offset,n);

    ending(result);

    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);

    //----------------------------------------------------------------------------

    //  Velocities
    iin=sizeof(float)*(p->pointnum2D)*3;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol1a(b->P));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->sl_ipol2a(b->Q));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);

	ffn=float(p->sl_ipol4(b->ws));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}


	//  Pressure
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->press));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  eddyv
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->eddyv));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  turbulence
    pturb->print_2D(p,b,pgc,buffer,m);

    //  Eta
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4eta_wd(p->wet,b->eta,b->bed));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

	//  Waterlevel
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->hp));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  wetdry
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=0.25*float((p->wet[IJ]+p->wet[Ip1J]+p->wet[IJp1]+p->wet[Ip1Jp1]));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    
    //  Breaking
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->breaking_print));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}

    //  test
    if(p->P23==1)
    {
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->test));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    //  fb
    if(p->P28==1)
    {
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->fs));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }
    
    //  Hs
    if(p->P110==1)
    {
	iin=sizeof(float)*(p->pointnum2D);
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->Hs));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
    }

    //  XYZ
	iin=sizeof(float)*(p->pointnum2D)*3;
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

        switch (p->P73)
        {
        case 0: default:
            ffn=p->sl_ipol4eta(p->wet,b->eta,b->bed)+p->wd;
            break;
        case 1:
            ffn=0.5*(b->hx(i,j)+b->hx(i,j+1)) + p->sl_ipol4(b->bed);
            break;
        case 2:
            ffn=0.5*(b->hy(i,j)+b->hy(i+1,j)) + p->sl_ipol4(b->bed);
            break;
        }
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
	}

    //  Connectivity
    iin=sizeof(int)*(p->polygon_sum)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(b->nodeval(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(b->nodeval(i,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(b->nodeval(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);


	// Triangle 2
	iin=int(b->nodeval(i-1,j-1))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(b->nodeval(i,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	iin=int(b->nodeval(i-1,j))-1;
	std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	}


    // Offset of Connectivity
    iin=sizeof(int)*(p->polygon_sum);
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

    auto stop = chrono::high_resolution_clock::now();
    auto duration = pgc->globalmax(chrono::duration_cast<chrono::microseconds>(stop - start).count());
    if(p->mpirank==0)
    std::cerr<<duration<<std::endl;

	++p->printcount;

}
