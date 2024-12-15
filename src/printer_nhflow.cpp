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

#include"printer_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sediment.h"
#include"nhflow_print_wsf.h"
#include"nhflow_vtp_fsf.h"
#include"nhflow_vtp_bed.h"
#include"nhflow_state.h"
#include"nhflow_print_wsf_theory.h"
#include"nhflow_print_wsfline.h"
#include"nhflow_print_wsfline_y.h"
#include"nhflow_print_runup_gage_x.h"
#include"nhflow_print_runup_max_gage_x.h"
#include"nhflow_vel_probe.h"
#include"nhflow_vel_probe_theory.h"
#include"nhflow_print_Hs.h"
#include"nhflow_turbulence.h"
#include"nhflow_force.h"
#include"nhflow_force_ale.h"
#include"bedshear_probe.h"
#include"bedshear_max.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<sstream>
#include<cstdio>
#include<cstring>

printer_nhflow::printer_nhflow(lexer* p, fdm_nhf *d, ghostcell *pgc)
{	
    switch (p->P10)
    {
        case 0: case 2:
            outputFormat = new vtk3D();
            break;
        case 1: default:
            outputFormat = new vtu3D();
            break;
        case 3:
            outputFormat = new vts3D();
            break;
    }

    if(p->I40==0)
    {
	p->printtime=0.0;
	p->sedprinttime=0.0;
	p->fsfprinttime=0.0;
	p->probeprinttime=0.0;
	p->stateprinttime=0.0;
    p->exportprinttime=0.0;
    }

	p->Darray(printtime_wT,p->P35);
    p->Iarray(printfsfiter_wI,p->P184);
    p->Darray(printfsftime_wT,p->P185);
    
    
    p->Iarray(printfsfiter_wI,p->P184);

	for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];

    for(int qn=0; qn<p->P185; ++qn)
	printfsftime_wT[qn]=p->P185_ts[qn];

    for(int qn=0; qn<p->P184; ++qn)
	printfsfiter_wI[qn]=p->P184_its[qn];


	printcount=0;

	// Create Folder
	if(p->mpirank==0)
	outputFormat->folder("NHFLOW");
    
    pwsf=new nhflow_print_wsf(p,d);

    pwsf_theory=new nhflow_print_wsf_theory(p,d,pgc);

    pwsfline=new nhflow_print_wsfline(p,d,pgc);

    pwsfline_y=new nhflow_print_wsfline_y(p,d,pgc);
    
    if(p->P65>0)
    pvel=new nhflow_vel_probe(p,d);
    
    if(p->P66>0)
    pveltheo=new nhflow_vel_probe_theory(p,d);
    
    prunupx=new nhflow_print_runup_gage_x(p,d,pgc);
    
    prunupmaxx=new nhflow_print_runup_max_gage_x(p,d,pgc);
    
    if(p->P40>0)
	pstate=new nhflow_state(p,d,pgc);
    
    if(p->P81>0)
	pforce = new nhflow_force*[p->P81];
    
    for(n=0;n<p->P81;++n)
	pforce[n]=new nhflow_force(p,d,pgc,n);
    
    if(p->P85>0)
    {
	pforce_ale = new nhflow_force_ale*[p->P85];
    
    for(n=0;n<p->P85;++n)
	pforce_ale[n]=new nhflow_force_ale(p,d,pgc,n);
    
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->ini(p,d,pgc);
    }
    
    if(p->P110==1)
    phs = new nhflow_print_Hs(p,d->Hs);
    
    if(p->P180==1)
	pfsf = new nhflow_vtp_fsf(p,d,pgc);

    pbed = new nhflow_vtp_bed(p,d,pgc);

    if(p->P125>0)
    pbedshear = new bedshear_probe(p,pgc);

    if(p->P126==1)
    pbedshearmax = new bedshear_max(p,pgc);
    
}

printer_nhflow::~printer_nhflow()
{
}

void printer_nhflow::start(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow *pflow, nhflow_turbulence *pnhfturb, sediment *psed)
{
    // Gages
	if(p->P51>0)
	pwsf->height_gauge(p,d,pgc,d->eta);

    if(p->P50>0)
    pwsf_theory->height_gauge(p,d,pgc,pflow);
    
    if(p->P110==1)
    phs->start(p,pgc,d->eta,d->Hs);
    
    if(p->P133>0)
	prunupx->start(p,d,pgc,pflow,d->eta);
    
    if(p->P134>0)
	prunupmaxx->start(p,d,pgc,pflow,d->eta);
    
    if(p->P65>0)
	pvel->start(p,d,pgc);
    
    if(p->P66>0)
	pveltheo->start(p,d,pgc,pflow);
    
    pfsf->preproc(p,d,pgc);

    // Print out based on iteration
    if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P20>0)
    {
    print_vtu(p,d,pgc,pnhfturb,psed);
    }

    // Print out based on time
    if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0) || (p->count==0 &&  p->P30>0.0))
    {
    print_vtu(p,d,pgc,pnhfturb,psed);
    
    p->printtime+=p->P30;
    }

    // Print out based on time interval
    if(p->P35>0)
    for(int qn=0; qn<p->P35; ++qn)
    if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
    {
    print_vtu(p,d,pgc,pnhfturb,psed);

    printtime_wT[qn]+=p->P35_dt[qn];
    }

    // Print FSF
    if(((p->count%p->P181==0 && p->P182<0.0 && p->P180==1 )|| (p->count==0 &&  p->P182<0.0 && p->P180==1)) && p->P181>0)
    {
    pfsf->start(p,d,pgc,psed);
    
    if(p->S10>0)
    pbed->start(p,d,pgc,psed);
    }


    if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
    {
    pfsf->start(p,d,pgc,psed);
    
    if(p->S10>0)
    pbed->start(p,d,pgc,psed);
    
    p->fsfprinttime+=p->P182;
    }

    if(p->P180==1 && p->P184>0)
    for(int qn=0; qn<p->P184; ++qn)
    if(p->count%p->P184_dit[qn]==0 && p->count>=p->P184_its[qn] && p->count<=(p->P184_ite[qn]))
    {
    pfsf->start(p,d,pgc,psed);
    
    if(p->S10>0)
    pbed->start(p,d,pgc,psed);
    }

    if(p->P180==1 && p->P185>0)
    for(int qn=0; qn<p->P185; ++qn)
    if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P185_ts[qn] && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
    {
    pfsf->start(p,d,pgc,psed);
    
    if(p->S10>0)
    pbed->start(p,d,pgc,psed);

    printfsftime_wT[qn]+=p->P185_dt[qn];
    }

    // Print BED
    if(p->count==0 && p->S10==0)
    pbed->start(p,d,pgc,psed);


    // Gages
    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,d,pgc,pflow,d->eta);

    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,d,pgc,pflow,d->eta);


    // Print state out based on iteration
    if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && (p->P46==0 || (p->count>=p->P46_is && p->count<p->P46_ie)))
    {
    pstate->write(p,d,pgc);
    }

    // Print sate out based on time
    if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0 && (p->P47==0 || (p->count>=p->P47_ts && p->count<<p->P47_te)))
    {
    pstate->write(p,d,pgc);

    p->stateprinttime+=p->P42;
    }
    
    // Force box
    if((p->count==0 || p->count==p->count_statestart) && p->P81>0)
	for(n=0;n<p->P81;++n)
	pforce[n]->ini(p,d,pgc);

	if(p->count>1 && p->P81>0)
	for(n=0;n<p->P81;++n)
	pforce[n]->start(p,d,pgc);
    
    
    // ALE force    
    if(p->count>0)
    if(p->count%p->P80==0)
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->start(p,d,pgc);

    /*
    if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
    p->probeprinttime+=p->P55;

    if(p->P59==1)
    pbreaklog->write(p,d,pgc);
    */

    // sediment probes
    if(p->P125>0)
	pbedshear->bedshear_gauge(p,pgc,psed);

    if(p->P126==1)
    pbedshearmax->bedshear_maxval(p,pgc,psed);
}

void printer_nhflow::print_stop(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow *pflow, nhflow_turbulence *pnhfturb, sediment *psed)
{
    print_vtu(p,d,pgc,pnhfturb,psed);
    
    if(p->S10>0)
    pbed->start(p,d,pgc,psed);
    
    if(p->P180==1)
    pfsf->start(p,d,pgc,psed);
    
}

void printer_nhflow::print_vtu(lexer* p, fdm_nhf *d, ghostcell* pgc, nhflow_turbulence *pnhfturb, sediment *psed)
{
    if(p->P10==1||p->P10==3)
    {
        /*
        - U, V, W
        - P
        - test
        - breaking
        */
        
        SLICELOOP4
        {
        if(d->breaking(i,j)==1)
        d->breaking_print(i,j)=1.0;

        if(d->breaking(i,j)==0)
        d->breaking_print(i,j)=0.0;
        }
        
        //
        //pgc->gcsl_start4(p,d->WL,50);
        pgc->gcsl_start4(p,d->bed,50);
        pgc->gcsl_start4(p,d->breaking_print,50);
        pgc->start4V(p,d->test,50);
        //pgc->start4(p,d->test,1);
        
        pgc->dgcslpol(p,d->WL,p->dgcsl4,p->dgcsl4_count,14);
        pgc->dgcslpol(p,d->breaking_print,p->dgcsl4,p->dgcsl4_count,14);
        pgc->dgcslpol(p,d->bed,p->dgcsl4,p->dgcsl4_count,14);

        d->WL.ggcpol(p);
        d->breaking_print.ggcpol(p);

        i=-1;
        j=-1;
        if(i+p->origin_i==-1 && j+p->origin_j==-1 )
        d->WL(i,j) = d->WL(i+1,j+1);


        //----------

        outputFormat->extent(p,pgc);
        if(p->mpirank==0)
        parallel(p,d,pgc,pnhfturb,psed);

        int num=0;
        if(p->P15==1)
        num = printcount;
        if(p->P15==2)
        num = p->count;
        outputFormat->fileName(name,"NHFLOW",num,p->mpirank+1);

        if(initial_print)
        {
            initial_print = false;

            n=0;

            offset[n]=0;
            ++n;

            // velocity
            offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
            ++n;

            // scalars

            // P
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            
            // k and eps
            pnhfturb->offset_vtu(p,offset,n);
            
            // omega_sig
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;

            // elevation
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            
            // test
            if(p->P23==1)
            {
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            }
            // Hs
            if(p->P110==1)
            {
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            }
            
            // solid
            if(p->P25==1)
            {
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            }
            
            if(p->P25==1 || p->P28==1)  
            {
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            }
            
            // floating
            if(p->P28==1)
            {
            offset[n]=offset[n-1]+4*(p->pointnum)+4;
            ++n;
            }

            // Format specific structure
            outputFormat->offset(p,offset,n);
        }
        //---------------------------------------------

        std::stringstream result;
        outputFormat->beginning(p,result);

        n=0;
        result<<"<PointData>\n";
        result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        result<<"<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
        
        pnhfturb->name_vtu(p,d,pgc,result,offset,n);
        
        result<<"<DataArray type=\"Float32\" Name=\"omega_sig\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
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


        if(p->P25==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"solid\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
        }
        
        if(p->P25==1 || p->P28==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"Heaviside\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
        }
        
        if(p->P28==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"floating\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
        }
        
        result<<"</PointData>\n";

        outputFormat->ending(result,offset,n);
        result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

        m=result.str().length();
        buffer.resize(m+offset[n]+27);
        std::memcpy(&buffer[0],result.str().data(),m);

    //----------------------------------------------------------------------------
        

    //  Velocities
        iin=3*4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->U[IJK]+d->U[Ip1JK]+d->U[IJKp1]+d->U[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->U[IJK]+d->U[Ip1JK]+d->U[IJp1K]+d->U[Ip1Jp1K]
                    +  d->U[IJKp1]+d->U[Ip1JKp1]+d->U[IJp1Kp1]+d->U[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);


        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->V[IJK]+d->V[Ip1JK]+d->V[IJKp1]+d->V[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->V[IJK]+d->V[Ip1JK]+d->V[IJp1K]+d->V[Ip1Jp1K]
                    +  d->V[IJKp1]+d->V[Ip1JKp1]+d->V[IJp1Kp1]+d->V[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);


        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->W[IJK]+d->W[Ip1JK]+d->W[IJKp1]+d->W[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->W[IJK]+d->W[Ip1JK]+d->W[IJp1K]+d->W[Ip1Jp1K]
                    +  d->W[IJKp1]+d->W[Ip1JKp1]+d->W[IJp1Kp1]+d->W[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }

    //  P
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.5*(d->P[FIJKp1]+d->P[FIm1JKp1]));
        j=jj;
        }
            
        if(p->j_dir==1)
        ffn=float(0.25*(d->P[FIJKp1]+d->P[FIm1JKp1] + d->P[FIJm1Kp1]+d->P[FIm1Jm1Kp1]));

        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
    
    //  kin and eps
        pnhfturb->print_3D(p,d,pgc,buffer,m);
        
    //  Omega_sig
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float((d->omegaF[FIJKp1]));
        j=jj;
        }

        if(p->j_dir==1)
        ffn=float(0.5*(d->omegaF[FIJKp1]+d->omegaF[FIJp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }

    //  elevation
        iin=4*(p->pointnum)*3;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        ffn=float(p->ZN[KP1]*d->WL(i,j) + d->bed(i,j));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }

    //  test
        if(p->P23==1)
        {
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->test[IJK]+d->test[Ip1JK]+d->test[IJKp1]+d->test[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->test[IJK]+d->test[Ip1JK]+d->test[IJp1K]+d->test[Ip1Jp1K]
                    +  d->test[IJKp1]+d->test[Ip1JKp1]+d->test[IJp1Kp1]+d->test[Ip1Jp1Kp1]));
                    
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
        }
        
    //  Hs
        if(p->P110==1)
        {
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        ffn=float(p->sl_ipol4(d->Hs));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
        }
        
    //  solid
        if(p->P25==1)
        {
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->SOLID[IJK]+d->SOLID[Ip1JK]+d->SOLID[IJKp1]+d->SOLID[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->SOLID[IJK]+d->SOLID[Ip1JK]+d->SOLID[IJp1K]+d->SOLID[Ip1Jp1K]
                    +  d->SOLID[IJKp1]+d->SOLID[Ip1JKp1]+d->SOLID[IJp1Kp1]+d->SOLID[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
        }
        
        if(p->P25==1 || p->P28==1)
        {
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->FHB[IJK]+d->FHB[Ip1JK]+d->FHB[IJKp1]+d->FHB[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->FHB[IJK]+d->FHB[Ip1JK]+d->FHB[IJp1K]+d->FHB[Ip1Jp1K]
                    +  d->FHB[IJKp1]+d->FHB[Ip1JKp1]+d->FHB[IJp1Kp1]+d->FHB[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
        }
        
    //  floating
        if(p->P28==1)
        {
        iin=4*(p->pointnum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        TPLOOP
        {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->FB[IJK]+d->FB[Ip1JK]+d->FB[IJKp1]+d->FB[Ip1JKp1]));
        j=jj;
        }
        
        if(p->j_dir==1)
        ffn=float(0.125*(d->FB[IJK]+d->FB[Ip1JK]+d->FB[IJp1K]+d->FB[Ip1Jp1K]
                    +  d->FB[IJKp1]+d->FB[Ip1JKp1]+d->FB[IJp1Kp1]+d->FB[Ip1Jp1Kp1]));
        
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
        }
        }

        // -----------------------

        outputFormat->structureWrite(p,d,buffer,m);

        // Open File
        FILE* file = fopen(name, "w");
        fwrite(buffer.data(), buffer.size(), 1, file);
        fclose(file);

        ++printcount;
    }
}
