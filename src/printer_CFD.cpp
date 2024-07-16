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

#include"printer_CFD.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"solver.h"
#include"print_wsf.h"
#include"print_wsf_theory.h"
#include"print_wsfline_x.h"
#include"print_wsfline_y.h"
#include"force.h"
#include"vorticity_f.h"
#include"vorticity_void.h"
#include"probe_point.h"
#include"probe_pressure.h"
#include"probe_line.h"
#include"bedprobe_point.h"
#include"bedprobe_max.h"
#include"ioflow.h"
#include"data.h"
#include"concentration.h"
#include"gage_discharge_x.h"
#include"gage_discharge_window_x.h"
#include"fsf_vtp.h"
#include"topo_vtp.h"
#include"cfd_state.h"
#include"bedshear_probe.h"
#include"bedshear_max.h"
#include"bedprobe_line_x.h"
#include"bedprobe_line_y.h"
#include"probe_vel.h"
#include"probe_vel_theory.h"
#include"multiphase.h"
#include"sediment.h"
#include"sloshing_force.h"
#include"print_porous.h"
#include"flowfile_out.h"
#include"print_averaging_f.h"
#include"print_averaging_v.h"
#include<sys/stat.h>
#include<sys/types.h>

#include "vtks.h"

printer_CFD::printer_CFD(lexer* p, fdm *a, ghostcell *pgc) : nodefill(p), eta(p)
{
    switch (p->P10)
    {
        case 0:
            outputFormat = new vtk3D();
            break;
        case 1: default:
            outputFormat = new vtu3D();
            break;
        case 2:
            outputFormat = new vtr3D();
            break;
        case 3:
            outputFormat = new vts3D();
            break;
    }

    // if(p->F50==1)
	// gcval_phi=51;

	// if(p->F50==2)
	// gcval_phi=52;

	// if(p->F50==3)
	// gcval_phi=53;

	// if(p->F50==4)
	// gcval_phi=54;

	// if(p->F50==1)
	// gcval_phiext=61;

	// if(p->F50==2)
	// gcval_phiext=62;

	// if(p->F50==3)
	// gcval_phiext=63;

	// if(p->F50==4)
	// gcval_phiext=64;

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
    // p->Iarray(printfsfiter_wI,p->P184);
    p->Darray(printfsftime_wT,p->P185);

	for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];

    for(int qn=0; qn<p->P185; ++qn)
	printfsftime_wT[qn]=p->P185_ts[qn];

    // for(int qn=0; qn<p->P184; ++qn)
	// printfsfiter_wI[qn]=p->P184_its[qn];

	pwsf=new print_wsf(p,a,pgc,0);
	pwsf_theory=new print_wsf_theory(p,a,pgc,0);
	pwsfline_x=new print_wsfline_x(p,a,pgc);
	pwsfline_y=new print_wsfline_y(p,a,pgc);
	pprobe = new probe_point(p,a,pgc);
    ppressprobe = new probe_pressure(p,a,pgc);
	pline = new probe_line(p,a,pgc);
	pq = new gage_discharge_x(p,a,pgc);
    pqw = new gage_discharge_window_x(p,a,pgc);
    
    if(p->P21==0)
    pmean = new print_averaging_v(p,a,pgc);
    
    if(p->P21==1)
    pmean = new print_averaging_f(p,a,pgc);

	if(p->P180==1)
	pfsf = new fsf_vtp(p,a,pgc);
    
    if(p->P190==1)
	ptopo = new topo_vtp(p,a,pgc);
    
    if(p->P65>0)
    pvel=new probe_vel(p,a);
    
    if(p->P66>0)
    pveltheo=new probe_vel_theory(p,a);

	if(p->P75==0)
	pvort = new vorticity_void(p,a);

	if(p->P75==1)
	pvort = new vorticity_f(p,a);

    if(p->P81>0)
    {
        P81 = p->P81;
	    pforce = new force*[P81];
    }

	if(p->P121>0)
	pbedpt = new bedprobe_point(p,a,pgc);

	if(p->P122>0)
	pbedmax = new bedprobe_max(p,a,pgc);

	if(p->P123>0)
	pbedlinex=new bedprobe_line_x(p,a,pgc);

	if(p->P124>0)
	pbedliney=new bedprobe_line_y(p,a,pgc);

	if(p->P125>0)
	pbedshear = new bedshear_probe(p,a,pgc);

	if(p->P126>0)
	pbedshearmax = new bedshear_max(p,a,pgc);

    for(n=0;n<p->P81;++n)
	pforce[n]=new force(p,a,pgc,n);

	if(p->P40>0)
	pstate=new cfd_state(p,a,pgc);

    if(p->P101>0)
	pslosh=new sloshing_force(p,a,pgc);

	if(p->B270>0 || p->B274>0 || p->B281>0 || p->B282>0 || p->B291>0 || p->B310>0 || p->B321>0 || p->B322>0 || p->B311>0)
	{
	ppor=new print_porous(p,a,pgc);
	ppor->start(p,a,pgc);
	}

    if(p->P230>0)
    pflowfile = new flowfile_out(p,a,pgc);


	p->printcount=0;

	// Create Folder
	if(p->mpirank==0)
    {
        outputFormat->folder("CFD");
    }
}

printer_CFD::~printer_CFD()
{
    delete pwsf;
	delete pwsf_theory;
    delete pwsfline_x;
	delete pwsfline_y;
    for (int n=0; n<P81; ++n)
    {
        delete pforce[n];
    }
    delete pforce;
    delete pforce;
    delete pvort;
	delete pprobe;
    delete ppressprobe;
	delete pline;
	delete pbedpt;
	delete pbedlinex;
	delete pbedliney;
	delete pbedmax;
	delete pbedshear;
	delete pbedshearmax;
	delete pq;
    delete pqw;
	delete pfsf;
    delete ptopo;
	delete pstate;
    delete pslosh;
	delete ppor;
    delete pexport;
    delete pflowfile;
    delete pmean;
    delete pvel;
    delete pveltheo;

    delete outputFormat;
}

void printer_CFD::start(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
	pgc->gcparax4a(p,a->phi,5);
	
	pmean->averaging(p,a,pgc,pheat);

	// Print out based on iteration
	if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P20>0)
	{
	print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
	}

	// Print out based on time
	if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0) || (p->count==0 &&  p->P30>0.0))
	{
	print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

	p->printtime+=p->P30;
	}

	// Print out based on sediment time
	if((p->sedtime>p->sedprinttime && p->P34>0.0 && p->P30<0.0) || (p->count==0 &&  p->P34>0.0))
	{
	print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

	p->sedprinttime+=p->P34;
	}

	// Print out based on time interval
	if(p->P10!=0 && p->P35>0)
	for(int qn=0; qn<p->P35; ++qn)
	if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
	{
	print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

	printtime_wT[qn]+=p->P35_dt[qn];
	}


	if((p->P62>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P62>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
	pline->start(p,a,pgc,pturb);


	if(p->P50>0)
	pwsf_theory->height_gauge(p,a,pgc,pflow,a->phi);

	if(p->P51>0)
	pwsf->height_gauge(p,a,pgc,a->phi);

	if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
	pwsfline_x->wsfline(p,a,pgc,pflow);

	if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
	pwsfline_y->wsfline(p,a,pgc,pflow);

	if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
	{
	p->probeprinttime+=p->P55;
	}

	if(p->P61>0)
	pprobe->start(p,a,pgc,pturb);
	
	if(p->P64>0)
	ppressprobe->start(p,a,pgc,pturb);
	
	if(p->P65>0)
	pvel->start(p,a,pgc);
	
	if(p->P66>0)
	pveltheo->start(p,a,pgc,pflow);

	if(p->P167>0)
	pq->start(p,a,pgc);
	
	if(p->P168>0)
	pqw->start(p,a,pgc);

	if((p->count==0 || p->count==p->count_statestart) && p->P81>0)
	for(n=0;n<p->P81;++n)
	pforce[n]->ini(p,a,pgc);

	if(p->count>1 && p->P81>0)
	for(n=0;n<p->P81;++n)
	pforce[n]->start(p,a,pgc);

	if(p->P101>0)
	pslosh->start(p,a,pgc);

	// sediment probes
	if(((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47) ) && p->S10>0)
	if((p->S42==1 && p->count%p->S44==0 && p->sediter%p->P120==0) || (p->S42==2 && p->simtime>=p->sedsimtime && p->sediter%p->P120==0) || (p->S42==3  && p->simtime/p->wT>=p->sedwavetime && p->sediter%p->P120==0))
	{
	if(p->P121>0)
	pbedpt->bed_gauge(p,a,pgc);

	if(p->P122>0)
	pbedmax->bed_max(p,a,pgc);

	if(p->P123>0)
	pbedlinex->start(p,a,pgc,pflow);

	if(p->P124>0)
	pbedliney->start(p,a,pgc,pflow);

	if(p->P125>0)
	pbedshear->bedshear_gauge(p,a,pgc,psed);

	if(p->P126>0)
	pbedshearmax->bedshear_maxval(p,a,pgc,psed);
	}

	// Multiphase
	pmp->print_file(p,a,pgc);

	// Print FSF
	if(((p->count%p->P181==0 && p->P182<0.0 && p->P180==1 )|| (p->count==0 &&  p->P182<0.0 && p->P180==1)) && p->P181>0)
	pfsf->start(p,a,pgc);

	if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
	{
	pfsf->start(p,a,pgc);
	p->fsfprinttime+=p->P182;
	}

	if(p->P180==1 && p->P184>0)
	for(int qn=0; qn<p->P184; ++qn)
	if(p->count%p->P184_dit[qn]==0 && p->count>=p->P184_its[qn] && p->count<=(p->P184_ite[qn]))
	{
	pfsf->start(p,a,pgc);
	}

	if(p->P180==1 && p->P185>0)
	for(int qn=0; qn<p->P185; ++qn)
	if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P185_ts[qn] && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
	{
	pfsf->start(p,a,pgc);

	printfsftime_wT[qn]+=p->P185_dt[qn];
	}
	
	// Print TOPO
	if(((p->count%p->P191==0 && p->P182<0.0 && p->P190==1 )|| (p->count==0 &&  p->P192<0.0 && p->P190==1)) && p->P191>0)
	ptopo->start(p,a,pgc,psed);

	if((p->simtime>p->fsfprinttime && p->P192>0.0 && p->P190==1) || (p->count==0 &&  p->P192>0.0))
	{
	ptopo->start(p,a,pgc,psed);
	p->fsfprinttime+=p->P192;
	}

	if(p->P190==1 && p->P194>0)
	for(int qn=0; qn<p->P194; ++qn)
	if(p->count%p->P194_dit[qn]==0 && p->count>=p->P194_its[qn] && p->count<=(p->P194_ite[qn]))
	{
	ptopo->start(p,a,pgc,psed);
	}

	if(p->P190==1 && p->P195>0)
	for(int qn=0; qn<p->P195; ++qn)
	if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P195_ts[qn] && p->simtime<=(p->P195_te[qn]+0.5*p->P195_dt[qn]))
	{
	ptopo->start(p,a,pgc,psed);

	printfsftime_wT[qn]+=p->P195_dt[qn];
	}

	if(p->P230>0)
	pflowfile->start(p,a,pgc,pturb);

	// Print state out based on iteration
	if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && p->P41>0 && (p->P46==0 || (p->count>=p->P46_is && p->count<<p->P46_ie)))
	{
	pstate->write(p,a,pgc,pturb,psed);
	}

	// Print state out based on time
	if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0 && (p->P47==0 || (p->count>=p->P47_ts && p->count<<p->P47_te)))
	{
	pstate->write(p,a,pgc,pturb,psed);

	p->stateprinttime+=p->P42;
	}

}

void printer_CFD::print_stop(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    print_vtk(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);
    
}

void printer_CFD::print_vtk(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P180==1)
	pfsf->start(p,a,pgc);
    
    print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
}

void printer_CFD::print3D(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P10!=0)
    {
        pgc->start4a(p,a->test,1);
        pgc->start1(p,a->u,110);
        pgc->start2(p,a->v,111);
        pgc->start3(p,a->w,112);


        pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
        pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
        pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
        pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
        pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
        pgc->dgcpol4(p,a->phi,14);
        pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
        pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
        pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
        //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

        a->u.ggcpol(p);
        a->v.ggcpol(p);
        a->w.ggcpol(p);
        a->press.ggcpol(p);
        a->eddyv.ggcpol(p);
        a->phi.ggcpol(p);
        a->conc.ggcpol(p);
        a->ro.ggcpol(p);
        a->visc.ggcpol(p);
        a->phi.ggcpol(p);
        a->fb.ggcpol(p);
        a->fbh4.ggcpol(p);
        //a->test.ggcpol(p);
        

        pgc->gcparacox(p,a->phi,50);
        pgc->gcparacox(p,a->phi,50);

        pgc->gcparacox(p,a->topo,150);
        pgc->gcparacox(p,a->topo,150);
        
        //pgc->start4a(p,a->topo,159);

        pgc->gcperiodicx(p,a->press,4);

        double XN[p->gknox+2];
        double YN[p->gknoy+2];
        double ZN[p->gknoz+2];

        double eddyv[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        double phi[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        
        int recvcounts[p->mpi_size],displs[p->mpi_size];
        int recvcount = p->knox+1;
        pgc->gather_int(&recvcount,1,recvcounts,1);
        int disp = p->origin_i+1;
        pgc->gather_int(&disp,1,displs,1);
        pgc->gatherv_double(p->XN+marge,p->knox+1, XN,recvcounts,displs);

        recvcount = p->knoy+1;
        pgc->gather_int(&recvcount,1,recvcounts,1);
        disp = p->origin_j+1;
        pgc->gather_int(&disp,1,displs,1);
        pgc->gatherv_double(p->YN+marge,p->knoy+1, YN,recvcounts,displs);

        recvcount = p->knoz+1;
        pgc->gather_int(&recvcount,1,recvcounts,1);
        disp = p->origin_k+1;
        pgc->gather_int(&disp,1,displs,1);
        pgc->gatherv_double(p->ZN+marge,p->knoz+1, ZN,recvcounts,displs);


        int neibours[6];
        neibours[0]=p->nb1;
        neibours[1]=p->nb2;
        neibours[2]=p->nb3;
        neibours[3]=p->nb4;
        neibours[4]=p->nb5;
        neibours[5]=p->nb6;
        int* gneibours;
        if(p->mpirank==0)
            gneibours = (int *)malloc(p->mpi_size*6*sizeof(int));
        pgc->gather_int(neibours,6,gneibours,6);

        int iextent[6];
        iextent[0]=p->origin_i;
        iextent[1]=p->origin_i+p->knox;
        iextent[2]=p->origin_j;
        iextent[3]=p->origin_j+p->knoy;
        iextent[4]=p->origin_k;
        iextent[5]=p->origin_k+p->knoz;

        int* piextent;
        if ( p->mpirank == 0)
        piextent = (int *)malloc(p->mpi_size*6*sizeof(int));
        pgc->gather_int(iextent,6,piextent,6);

        int localSendCount=0;
        localSendCount=(p->knox+2*p->margin)*(p->knoy+2*p->margin)*(p->knoz+2*p->margin);
        int globalSendCounts[p->mpi_size];
        pgc->gather_int(&localSendCount,1,globalSendCounts,1);

        double* pressGlobal;
        double* uvelGlobal;
        double* vvelGlobal;
        double* wvelGlobal;
        double* topoGlobal;
        double* phiGlobal;
        int counter;
        if(p->mpirank==0)
        {
            counter = globalSendCounts[0];
            displs[0]=0;
            for(int i=1;i<p->mpi_size;++i)
            {
                displs[i] = displs[i-1] + globalSendCounts[i-1];
                counter += globalSendCounts[i];
            }
            pressGlobal = (double *)malloc(counter*sizeof(double));
            uvelGlobal = (double *)malloc(counter*sizeof(double));
            vvelGlobal = (double *)malloc(counter*sizeof(double));
            wvelGlobal = (double *)malloc(counter*sizeof(double));
            topoGlobal = (double *)malloc(counter*sizeof(double));
            phiGlobal = (double *)malloc(counter*sizeof(double));
        }
        pgc->gatherv_double(a->press.V,localSendCount,pressGlobal,globalSendCounts,displs);
        pgc->gatherv_double(a->u.V,localSendCount,uvelGlobal,globalSendCounts,displs);
        pgc->gatherv_double(a->v.V,localSendCount,vvelGlobal,globalSendCounts,displs);
        pgc->gatherv_double(a->w.V,localSendCount,wvelGlobal,globalSendCounts,displs);
        pgc->gatherv_double(a->topo.V,localSendCount,topoGlobal,globalSendCounts,displs);
        pgc->gatherv_double(a->phi.V,localSendCount,phiGlobal,globalSendCounts,displs);

        if(p->mpirank==0)
        {
            double* press[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            double* uvel[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            double* vvel[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            double* wvel[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            double* topo[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            double* phi[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
            int indexL,indexLG;
            int kbegin,kend;
            int jbegin,jend;
            int ibegin,iend;
            for(int n=0;n<p->mpi_size;n++)
            {
                kbegin=-1;
                if(gneibours[4+6*n]>-2)
                kbegin=0;
                kend=piextent[5+6*n]-piextent[4+6*n]+1;
                if(gneibours[5+6*n]>-2)
                kend=piextent[5+6*n]-piextent[4+6*n];

                jbegin=-1;
                if(gneibours[2+6*n]>-2)
                jbegin=0;
                jend=piextent[3+6*n]-piextent[2+6*n]+1;
                if(gneibours[1+6*n]>-2)
                jend=piextent[3+6*n]-piextent[2+6*n];
                
                ibegin=-1;
                if(gneibours[0+6*n]>-2)
                ibegin=0;
                iend=piextent[1+6*n]-piextent[0+6*n]+1;
                if(gneibours[3+6*n]>-2)
                iend=piextent[1+6*n]-piextent[0+6*n];

                for(int k=kbegin;k<kend;++k)
                {
                    for(int j=jbegin;j<jend;++j)
                    {
                        for(int i=ibegin;i<iend;++i)
                        {
                            indexL = (i+p->margin)*(piextent[3+6*n]-piextent[2+6*n]+2*p->margin)*(piextent[5+6*n]-piextent[4+6*n]+2*p->margin) + (j+p->margin)*(piextent[5+6*n]-piextent[4+6*n]+2*p->margin) + k+p->margin;
                            indexL += displs[n];
                            indexLG = (k+1+piextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+piextent[2+6*n])*(p->gknox+2)+(i+1+piextent[0+6*n]);
                            press[indexLG]=&pressGlobal[indexL];
                            uvel[indexLG]=&uvelGlobal[indexL];
                            vvel[indexLG]=&vvelGlobal[indexL];
                            wvel[indexLG]=&wvelGlobal[indexL];
                            topo[indexLG]=&topoGlobal[indexL];
                            phi[indexLG]=&phiGlobal[indexL];
                        }
                    }
                }
            }
            int indexG;
            i=j=k=0;
            for(k=-1; k<1; ++k)
            for(j=-1; j<1; ++j)
            for(i=-1; i<1; ++i)
            {
                indexG = (k+1)*(p->gknox+2)*(p->gknoy+2)+(j+1)*(p->gknox+2)+(i+1);
                cout<<"("<<i<<","<<j<<","<<k<<")"<<a->press(i,j,k)<<":"<<*press[indexG]<<endl;
            }

            i=j=k=-1;
            XN[0]=p->XN[IP];
            YN[0]=p->YN[JP];
            ZN[0]=p->ZN[KP];            

            int cellNum=(p->gknox)*(p->gknoy)*(p->gknoz);
            int pointNum=(p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);

            int testOffset[300];
            int m=0;
            testOffset[m]=0;
            ++m;

            testOffset[m]=testOffset[m-1]+4+3*4*(pointNum);
            ++m;
            testOffset[m]=testOffset[m-1]+4+4*(pointNum);
            ++m;
            testOffset[m]=testOffset[m-1]+4+4*(pointNum);
            ++m;
            testOffset[m]=testOffset[m-1]+4+4*(pointNum);
            ++m;
            testOffset[m]=testOffset[m-1]+4+4*(pointNum);
            ++m;
            testOffset[m]=testOffset[m-1]+4+3*4*(cellNum);
            ++m;

            //x
            testOffset[m]=testOffset[m-1]+4+4*(p->gknox+1);
            ++m;
            //y
            testOffset[m]=testOffset[m-1]+4+4*(p->gknoy+1); 
            ++m;
            //z
            testOffset[m]=testOffset[m-1]+4+4*(p->gknoz+1);

            m=0;
            ofstream testFile;
            testFile.open("testFile.vtr",ios::binary);
            if(testFile.is_open())
            {
                testFile<<"<?xml version=\"1.0\"?>\n"
                <<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                <<"<RectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n"
                <<"<Piece Extent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\">\n"
                <<"<PointData>\n"
                <<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"</PointData>\n"
                <<"<CellData>\n"
                <<"\t<DataArray type=\"Float32\" Name=\"Debug\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<testOffset[m]<<"\" />\n";
                ++m;
                testFile<<"</CellData>\n"
                <<"<Coordinates>\n"
                <<"\t<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<testOffset[m]<<"\"/>\n";
                m++;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<testOffset[m]<<"\"/>\n";
                m++;
                testFile<<"\t<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<testOffset[m]<<"\"/>\n";
                m++;
                testFile<<"</Coordinates>\n"
                <<"</Piece>\n"
                <<"</RectilinearGrid>\n"
                <<"<AppendedData encoding=\"raw\">\n_";
                //  Velocities
                iin=3*4*(pointNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=0; k<p->gknoz+1; ++k)
                for(j=0; j<p->gknoy+1; ++j)
                for(i=0; i<p->gknox+1; ++i)
                {
                    ffn=float(p->ipol4press(uvel));//u
                    testFile.write((char*)&ffn, sizeof (float));

                    ffn=float(p->ipol4press(vvel));//v
                    testFile.write((char*)&ffn, sizeof (float));

                    ffn=float(p->ipol4press(wvel));//w
                    testFile.write((char*)&ffn, sizeof (float));
                }
                //  Pressure
                iin=4*(pointNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol4press(press));
                    testFile.write((char*)&ffn, sizeof (float));
                }
                //  EddyV
                iin=4*(pointNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=0; k<p->gknoz+1; ++k)
                for(j=0; j<p->gknoy+1; ++j)
                for(i=0; i<p->gknox+1; ++i)
                {
                ffn=float(p->ipol4_b(eddyv));//EddyV
                testFile.write((char*)&ffn, sizeof (float));
                }
                //  Phi
                iin=4*(pointNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol4phi(topo,phi));
                    testFile.write((char*)&ffn, sizeof (float));
                }
                //  Elevation
                iin=4*(pointNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=0; k<p->gknoz+1; ++k)
                for(j=0; j<p->gknoy+1; ++j)
                for(i=0; i<p->gknox+1; ++i)
                {
                    ffn=float(ZN[k]+(ZN[k+1]-ZN[k]));
                    testFile.write((char*)&ffn, sizeof (float));
                }

                //  Debug
                iin=3*4*(cellNum);
                testFile.write((char*)&iin, sizeof (int));
                for(k=0; k<p->gknoz; ++k)
                for(j=0; j<p->gknoy; ++j)
                for(i=0; i<p->gknox; ++i)
                {
                    ffn=float(*uvel[(k+1)*(p->gknox+2)*(p->gknoy+2)+(j+1)*(p->gknox+2)+(i+1)]);
                    testFile.write((char*)&ffn, sizeof (float));
                    ffn=float(*vvel[(k+1)*(p->gknox+2)*(p->gknoy+2)+(j+1)*(p->gknox+2)+(i+1)]);
                    testFile.write((char*)&ffn, sizeof (float));
                    ffn=float(*wvel[(k+1)*(p->gknox+2)*(p->gknoy+2)+(j+1)*(p->gknox+2)+(i+1)]);
                    testFile.write((char*)&ffn, sizeof (float));
                }

                // x
                iin=4*(p->gknox+1);
                testFile.write((char*)&iin, sizeof (int));
                for(i=1; i<p->gknox+2; ++i)
                {
                    ffn=float(XN[i]);
                    testFile.write((char*)&ffn, sizeof (float));
                }
                // y
                iin=4*(p->gknoy+1);
                testFile.write((char*)&iin, sizeof (int));
                for(j=1; j<p->gknoy+2; ++j)
                {
                    ffn=float(YN[j]);
                    testFile.write((char*)&ffn, sizeof (float));
                }
                // zle
                iin=4*(p->gknoz+1);
                testFile.write((char*)&iin, sizeof (int));
                for(k=1; k<p->gknoz+2; ++k)
                {
                    ffn=float(ZN[k]);
                    testFile.write((char*)&ffn, sizeof (float));
                }

                testFile<<"\n</AppendedData>\n"
                <<"</VTKFile>"<<flush;

                testFile.close();
            }
        }

        outputFormat->extent(p,pgc);
        if(p->mpirank==0)
        parallelData(a,p,pgc,pturb,pheat,pdata,pconc,pmp,psed);

        int num=0;
        if(p->P15==1)
        num = p->printcount;
        if(p->P15==2)
        num = p->count;
        int rank = p->mpirank+1;
        outputFormat->fileName(name,"CFD",num,rank);

        // Open File
        ofstream result;
        result.open(name, ios::binary);
        if(result.is_open())
        {
            n=0;

            offset[n]=0;
            ++n;

            // velocity
            offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
            ++n;
            
            pmean->offset_vtk(p,a,pgc,result,offset,n);

            // scalars
                // pressure
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                // k and eps
                pturb->offset_vtk(p,a,pgc,result,offset,n);
                // eddyv
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                // phi
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                // T
                pheat->offset_vtk(p,a,pgc,result,offset,n);
                // Multiphase
                pmp->offset_vtk(p,a,pgc,result,offset,n);
                // vorticity
                pvort->offset_vtk(p,a,pgc,result,offset,n);
                // data
                pdata->offset_vtk(p,a,pgc,result,offset,n);
                // concentration
                pconc->offset_vtk(p,a,pgc,result,offset,n);
                // rho
                if(p->P24==1 && p->F300==0)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // viscosity
                if(p->P71==1)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // VOF
                if(p->P72==1)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // Fi
                if(p->A10==4)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // conc
                if(p->P26==1)
                { 
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // topo
                if(p->P27==1)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // sediment bedlaod
                if(p->P76==1)
                psed->offset_vtk_bedload(p,pgc,result,offset,n);

                // sediment parameters 1
                if(p->P77==1)
                psed->offset_vtk_parameter1(p,pgc,result,offset,n);

                // sediment parameters 2
                if(p->P78==1)
                psed->offset_vtk_parameter2(p,pgc,result,offset,n);

                // bed shear stress
                if(p->P79>=1)
                psed->offset_vtk_bedshear(p,pgc,result,offset,n);

                // test
                if(p->P23==1)
                {
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
                // elevation
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                // solid
                if(p->P25==1)
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
                // walldist
                if(p->P29==1)
                {   
                offset[n]=offset[n-1]+4*(p->pointnum)+4;
                ++n;
                }
            // end scalars
            outputFormat->offset(p,offset,n);
            //---------------------------------------------

            outputFormat->beginning(p,result);
            
            if(p->P16==1)
            {
            result<<"<FieldData>"<<endl;
            result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
            result<<"</DataArray>"<<endl;
            result<<"</FieldData>"<<endl;
            }

            n=0;
            result<<"<PointData>"<<endl;
            result<<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            
            pmean->name_vtk(p,a,pgc,result,offset,n);

            result<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;

            pturb->name_vtk(p,a,pgc,result,offset,n);

            result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;

            pheat->name_vtk(p,a,pgc,result,offset,n);
            
            pmp->name_vtk(p,a,pgc,result,offset,n);

            pvort->name_vtk(p,a,pgc,result,offset,n);

            pdata->name_vtk(p,a,pgc,result,offset,n);

            pconc->name_vtk(p,a,pgc,result,offset,n);

            if(p->P24==1 && p->F300==0)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->P71==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"viscosity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }
            
            if(p->P72==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"VOF\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->A10==4)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->P26==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"ST_conc\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->P27==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"topo\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }
            
            if(p->P76==1)
            psed->name_vtk_bedload(p,pgc,result,offset,n);
            
            if(p->P77==1)
            psed->name_vtk_parameter1(p,pgc,result,offset,n);

            if(p->P78==1)
            psed->name_vtk_parameter2(p,pgc,result,offset,n);

            if(p->P79>=1)
            psed->name_vtk_bedshear(p,pgc,result,offset,n);

            if(p->P23==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            result<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;

            if(p->P25==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->P28==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"floating\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }

            if(p->P29==1)
            {
            result<<"\t<DataArray type=\"Float32\" Name=\"walldist\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
            ++n;
            }
            result<<"</PointData>"<<endl;
            outputFormat->ending(result,offset,n);
            //----------------------------------------------------------------------------
            result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

            //  Velocities
            iin=3*4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol1(a->u));
            result.write((char*)&ffn, sizeof (float));

            ffn=float(p->ipol2(a->v));
            result.write((char*)&ffn, sizeof (float));

            ffn=float(p->ipol3(a->w));
            result.write((char*)&ffn, sizeof (float));
            }

            //  time average flow parameters
            pmean->print_3D(p,a,pgc,result);

            //  Pressure
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4press(a->press)-p->pressgage);
            result.write((char*)&ffn, sizeof (float));
            }

            //  turbulence
            pturb->print_3D(p,a,pgc,result);

            //  eddyv
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->eddyv));
            result.write((char*)&ffn, sizeof (float));
            }

            //  phi
            nodefill4(p,a,pgc,a->phi,eta);
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            if(p->P18==1)
            ffn=float(p->ipol4phi(a,a->phi));
            if(p->P18==2)
            ffn = float(eta(i,j,k));
            result.write((char*)&ffn, sizeof (float));
            }

            //  T
            pheat->print_3D(p,a,pgc,result);
            
            //  Multiphase
            pmp->print_3D(p,a,pgc,result);

            //  Vorticity
            pvort->print_3D(p,a,pgc,result);

            //  Data
            pdata->print_3D(p,a,pgc,result);

            //  Concentration
            pconc->print_3D(p,a,pgc,result);

            //  density
            if(p->P24==1 && p->F300==0)
            {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->ro));
            result.write((char*)&ffn, sizeof (float));
            }
            }
            
            //  viscosity
            if(p->P71==1)
            {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4(a->visc));
            result.write((char*)&ffn, sizeof (float));
            }
            }
            
            //  VOF
            if(p->P72==1)
            {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4(a->vof));
            result.write((char*)&ffn, sizeof (float));
            }
            }

            //  Fi
            if(p->A10==4)
            {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4press(a->Fi));
            result.write((char*)&ffn, sizeof (float));
            }

            }

            if(p->P26==1)
            {
            //  conc
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4(a->conc));
            result.write((char*)&ffn, sizeof (float));
            }
            }
            
            if(p->P27==1)
            {
            //  topo
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->topo));
            //ffn = float(-a->bed(i,j)+p->ZN[KP1]);
            result.write((char*)&ffn, sizeof (float));
            }
            }
            
            //  sediment bedload
            if(p->P76==1)
            psed->print_3D_bedload(p,pgc,result);
            
            //  sediment parameter 1
            if(p->P77==1)
            psed->print_3D_parameter1(p,pgc,result);

            //  sediment parameter 2
            if(p->P78==1)
            psed->print_3D_parameter2(p,pgc,result);

            //  bed shear stress
            if(p->P79>=1)
            psed->print_3D_bedshear(p,pgc,result);

            //  test
            if(p->P23==1)
            {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->test));
            result.write((char*)&ffn, sizeof (float));
            }
            }

            //  elevation
            iin=4*(p->pointnum)*3;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->pos_z()+0.5*p->DZN[KP]);
            result.write((char*)&ffn, sizeof (float));
            }

            if(p->P25==1)
            {
            //  solid
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->solid));
            result.write((char*)&ffn, sizeof (float));
            }
            }

            if(p->P28==1)
            {
            //  floating
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->fb));
            result.write((char*)&ffn, sizeof (float));
            }
            }

            if(p->P29==1)
            {
            //  walldist
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
            ffn=float(p->ipol4_a(a->walld));
            result.write((char*)&ffn, sizeof (float));
            }
            }

            // -----------------------
            
            outputFormat->structureWrite(p,a,result);

            result.close();

            ++p->printcount;
        }
        else
        cout<<"Failed to open output file."<<endl;


        pgc->start1(p,a->u,114);
        pgc->start2(p,a->v,115);
        pgc->start3(p,a->w,116);

        pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
        pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
        pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
        pgc->start4a(p,a->topo,150);

        a->u.ggcpol(p);
        a->v.ggcpol(p);
        a->w.ggcpol(p);
    }
}