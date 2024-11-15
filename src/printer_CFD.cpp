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
#include<chrono>

#include "vtks.h"
#include "printMethods.h"

printer_CFD::printer_CFD(lexer* p, fdm *a, ghostcell *pgc) : nodefill(p), eta(p)
{
    // switch (p->P10)
    // {
    //     case 0:
    //         outputFormat = new vtk3D();
    //         break;
    //     case 1: default:
    //         outputFormat = new vtu3D();
    //         break;
    //     case 2:
    //         outputFormat = new vtr3D();
    //         break;
    //     case 3:
    //         outputFormat = new vts3D();
    //         break;
    // }

    switch (1)
    {
        case 0: default:
            outputMethod = new printMethodSeparated(p);
            break;
        case 1: 
            outputMethod = new printMethodCompact(p);
            break;
        case 2:
            outputMethod = new printMethodCompactMPI(p);
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
	// if(p->mpirank==0)
    // {
    //     outputFormat->folder("CFD");
    // }

    // setupCompactPrint(p,a,pgc);
    // setupCompactMPIPrint(p,a,pgc);
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
	    p->probeprinttime+=p->P55;

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
            pfsf->start(p,a,pgc);

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
	    pstate->write(p,a,pgc,pturb,psed);

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

// void printer_CFD::setupCompactPrint(lexer *p, fdm *a, ghostcell * pgc)
// {
//     if(p->mpirank==0)
//     {
// 	    mkdir("./REEF3D_CFD_VTRC",0777);

//         XN = (double *)malloc((p->gknox+2)*sizeof(double));
//         YN = (double *)malloc((p->gknoy+2)*sizeof(double));
//         ZN = (double *)malloc((p->gknoz+2)*sizeof(double));

//         gneibours = (int *)malloc(p->mpi_size*6*sizeof(int));
//         gextent = (int *)malloc(p->mpi_size*6*sizeof(int));

//         globalSendCounts = (int *)malloc(p->mpi_size*sizeof(int));
//         recvcounts = (int *)malloc(p->mpi_size*sizeof(int));
//         displs = (int *)malloc(p->mpi_size*sizeof(int));

//         cellNum=(p->gknox)*(p->gknoy)*(p->gknoz);
//         pointNum=(p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);

//         flag = (int **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(int*));
//         flag5 = (int **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(int*));

//         // ---------------------------------------------------------
//         // Allocate memory for data to be printed
//         // ---------------------------------------------------------
        
//         uvel = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         vvel = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         wvel = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         press = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         eddyv = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         phi = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));
//         topo = (double **)malloc(((p->gknox+2)*(p->gknoy+2)*(p->gknoz+2))*sizeof(double*));

//         // ---------------------------------------------------------
//         // Pre-calulate offsets
//         // ---------------------------------------------------------

//         m=0;
//         compactOffset[m]=0;
//         ++m;

//         //velocities
//         compactOffset[m]=compactOffset[m-1]+4+3*4*(pointNum);
//         ++m;
//         //pressure
//         compactOffset[m]=compactOffset[m-1]+4+4*(pointNum);
//         ++m;
//         //eddyv
//         compactOffset[m]=compactOffset[m-1]+4+4*(pointNum);
//         ++m;
//         //phi
//         compactOffset[m]=compactOffset[m-1]+4+4*(pointNum);
//         ++m;
//         //elevation
//         compactOffset[m]=compactOffset[m-1]+4+4*(pointNum);
//         ++m;

//         //x
//         compactOffset[m]=compactOffset[m-1]+4+4*(p->gknox+1);
//         ++m;
//         //y
//         compactOffset[m]=compactOffset[m-1]+4+4*(p->gknoy+1); 
//         ++m;
//         //z
//         compactOffset[m]=compactOffset[m-1]+4+4*(p->gknoz+1);
//     }
//     else
//     {
//         XN = nullptr;
//         YN = nullptr;
//         ZN = nullptr;

//         gneibours = nullptr;
//         gextent = nullptr;

//         globalSendCounts = nullptr;
//         recvcounts = nullptr;
//         displs = nullptr;

//         press = nullptr;
//         uvel = nullptr;
//         vvel = nullptr;
//         wvel = nullptr;
//         topo = nullptr;
//         phi = nullptr;
//         eddyv = nullptr;
//         flag = nullptr;
//         flag5 = nullptr;

//         cellNum=0;
//         pointNum=0;
//     }

//     int recvcount = p->knox+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     int disp = p->origin_i+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->XN+marge,p->knox+1, XN,recvcounts,displs);

//     recvcount = p->knoy+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     disp = p->origin_j+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->YN+marge,p->knoy+1, YN,recvcounts,displs);

//     recvcount = p->knoz+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     disp = p->origin_k+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->ZN+marge,p->knoz+1, ZN,recvcounts,displs);

//     if(p->mpirank==0)
//     {
//         i=j=k=-1;
//         XN[0]=p->XN[IP];
//         YN[0]=p->YN[JP];
//         ZN[0]=p->ZN[KP]; 
//     }


//     int neibours[6];
//     neibours[0]=p->nb1;
//     neibours[1]=p->nb2;
//     neibours[2]=p->nb3;
//     neibours[3]=p->nb4;
//     neibours[4]=p->nb5;
//     neibours[5]=p->nb6;
//     pgc->gather_int(neibours,6,gneibours,6);

//     int extent[6];
//     extent[0]=p->origin_i;
//     extent[1]=p->origin_i+p->knox;
//     extent[2]=p->origin_j;
//     extent[3]=p->origin_j+p->knoy;
//     extent[4]=p->origin_k;
//     extent[5]=p->origin_k+p->knoz;
//     pgc->gather_int(extent,6,gextent,6);

//     if(p->mpirank==0)
//         localSendCount=0;
//     else
//         localSendCount=(p->knox+2*p->margin)*(p->knoy+2*p->margin)*(p->knoz+2*p->margin);
//     pgc->gather_int(&localSendCount,1,globalSendCounts,1);

//     if(p->mpirank==0)
//     {
//         int counter = globalSendCounts[0];
//         displs[0]=0;
//         for(int i=1;i<p->mpi_size;++i)
//         {
//             displs[i] = displs[i-1] + globalSendCounts[i-1];
//             counter += globalSendCounts[i];
//         }
//         pressGlobal = (double *)malloc(counter*sizeof(double));
//         uvelGlobal = (double *)malloc(counter*sizeof(double));
//         vvelGlobal = (double *)malloc(counter*sizeof(double));
//         wvelGlobal = (double *)malloc(counter*sizeof(double));
//         topoGlobal = (double *)malloc(counter*sizeof(double));
//         phiGlobal = (double *)malloc(counter*sizeof(double));
//         eddyvGlobal = (double *)malloc(counter*sizeof(double));
//         flagGlobal = (int *)malloc(counter*sizeof(int));
//         flag5Global = (int *)malloc(counter*sizeof(int));
//     }
//     else
//     {
//         pressGlobal = nullptr;
//         uvelGlobal = nullptr;
//         vvelGlobal = nullptr;
//         wvelGlobal = nullptr;
//         topoGlobal = nullptr;
//         phiGlobal = nullptr;
//         eddyvGlobal = nullptr;
//         flagGlobal = nullptr;
//         flag5Global = nullptr;
//     }

//     if(p->mpirank==0)
//     {
//         int indexL,indexLG;
//         int kbegin,kend;
//         int jbegin,jend;
//         int ibegin,iend;
//         for(int n=0;n<p->mpi_size;n++)
//         {
//             kbegin=-1;
//             if(gneibours[4+6*n]>-2)
//                 kbegin=0;
//             kend=gextent[5+6*n]-gextent[4+6*n]+1;
//             if(gneibours[5+6*n]>-2)
//                 kend=gextent[5+6*n]-gextent[4+6*n];

//             jbegin=-1;
//             if(gneibours[2+6*n]>-2)
//                 jbegin=0;
//             jend=gextent[3+6*n]-gextent[2+6*n]+1;
//             if(gneibours[1+6*n]>-2)
//                 jend=gextent[3+6*n]-gextent[2+6*n];
            
//             ibegin=-1;
//             if(gneibours[0+6*n]>-2)
//                 ibegin=0;
//             iend=gextent[1+6*n]-gextent[0+6*n]+1;
//             if(gneibours[3+6*n]>-2)
//                 iend=gextent[1+6*n]-gextent[0+6*n];

//             if(n!=0)
//                 for(int k=kbegin;k<kend;++k)
//                 {
//                     for(int j=jbegin;j<jend;++j)
//                     {
//                         for(int i=ibegin;i<iend;++i)
//                         {
//                             indexL = (i+p->margin)*(gextent[3+6*n]-gextent[2+6*n]+2*p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + (j+p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + k+p->margin;
//                             indexL += displs[n];
//                             indexLG = (k+1+gextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+gextent[2+6*n])*(p->gknox+2)+(i+1+gextent[0+6*n]);
//                             press[indexLG]=&pressGlobal[indexL];
//                             uvel[indexLG]=&uvelGlobal[indexL];
//                             vvel[indexLG]=&vvelGlobal[indexL];
//                             wvel[indexLG]=&wvelGlobal[indexL];
//                             topo[indexLG]=&topoGlobal[indexL];
//                             phi[indexLG]=&phiGlobal[indexL];
//                             eddyv[indexLG]=&eddyvGlobal[indexL];
//                             flag[indexLG]=&flagGlobal[indexL];
//                             flag5[indexLG]=&flag5Global[indexL];
//                         }
//                     }
//                 }
//             else
//                 for(int k=kbegin;k<kend;++k)
//                 {
//                     for(int j=jbegin;j<jend;++j)
//                     {
//                         for(int i=ibegin;i<iend;++i)
//                         {
//                             indexL = (i+p->margin)*(gextent[3+6*n]-gextent[2+6*n]+2*p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + (j+p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + k+p->margin;
//                             indexLG = (k+1+gextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+gextent[2+6*n])*(p->gknox+2)+(i+1+gextent[0+6*n]);
//                             press[indexLG]=&a->press.V[indexL];
//                             uvel[indexLG]=&a->u.V[indexL];
//                             vvel[indexLG]=&a->v.V[indexL];
//                             wvel[indexLG]=&a->w.V[indexL];
//                             topo[indexLG]=&a->topo.V[indexL];
//                             phi[indexLG]=&a->phi.V[indexL];
//                             eddyv[indexLG]=&a->eddyv.V[indexL];
//                             flag[indexLG]=&p->flag[indexL];
//                             flag5[indexLG]=&p->flag5[indexL];
//                         }
//                     }
//                 }
//         }
//     }
// }

// void printer_CFD::print3Dcompact(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
// {
//     if(p->P10==2)
//     {
//         pgc->gatherv_double(a->press.V,localSendCount,pressGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->u.V,localSendCount,uvelGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->v.V,localSendCount,vvelGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->w.V,localSendCount,wvelGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->topo.V,localSendCount,topoGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->phi.V,localSendCount,phiGlobal,globalSendCounts,displs);
//         pgc->gatherv_double(a->eddyv.V,localSendCount,eddyvGlobal,globalSendCounts,displs);
//         pgc->gatherv_int(p->flag,localSendCount,flagGlobal,globalSendCounts,displs);
//         pgc->gatherv_int(p->flag5,localSendCount,flag5Global,globalSendCounts,displs);

//         if(p->mpirank==0)
//         {
//             int num=0;
//             if(p->P15==1)
//                 num = p->printcount-1; // temporary fix
//             if(p->P15==2)
//                 num = p->count;
//             sprintf(name,"./REEF3D_CFD_VTRC/REEF3D-CFD-%08i.vtr",num);
//             m=0;
//             ofstream result;
//             result.open(name,ios::binary);
//             if(result.is_open())
//             {
//                 result<<"<?xml version=\"1.0\"?>\n"
//                 <<"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
//                 <<"<RectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
//                 if(p->P16==1)
//                 {
//                 result<<"<FieldData>"<<endl;
//                 result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<std::setprecision(7)<<p->simtime<<"</DataArray>"<<endl;
//                 result<<"</FieldData>"<<endl;
//                 }
//                 result<<"<Piece Extent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\">\n";
//                 result<<"<PointData>\n"
//                 <<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 ++m;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 ++m;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 ++m;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 ++m;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 ++m;
//                 result<<"</PointData>\n"
//                 // <<"<CellData>\n"
//                 // <<"\t<DataArray type=\"Float32\" Name=\"Debug\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<compactOffset[m]<<"\" />\n";
//                 // ++m;
//                 // result<<"</CellData>\n"
//                 <<"<Coordinates>\n"
//                 <<"\t<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<compactOffset[m]<<"\"/>\n";
//                 m++;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<compactOffset[m]<<"\"/>\n";
//                 m++;
//                 result<<"\t<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<compactOffset[m]<<"\"/>\n";
//                 m++;
//                 result<<"</Coordinates>\n"
//                 <<"</Piece>\n"
//                 <<"</RectilinearGrid>\n"
//                 <<"<AppendedData encoding=\"raw\">\n_";
//                 //  Velocities
//                 iin=3*4*(pointNum);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=-1; k<p->gknoz; ++k)
//                     for(j=-1; j<p->gknoy; ++j)
//                         for(i=-1; i<p->gknox; ++i)
//                         {
//                             ffn=float(p->ipol1(uvel,flag,flag5));//u
//                             result.write((char*)&ffn, sizeof (float));

//                             ffn=float(p->ipol2(vvel,flag,flag5));//v
//                             result.write((char*)&ffn, sizeof (float));

//                             ffn=float(p->ipol3(wvel,flag,flag5));//w
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                 //  Pressure
//                 iin=4*(pointNum);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=-1; k<p->gknoz; ++k)
//                     for(j=-1; j<p->gknoy; ++j)
//                         for(i=-1; i<p->gknox; ++i)
//                         {
//                             ffn=float(p->ipol4press(press));
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                 //  EddyV
//                 iin=4*(pointNum);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=-1; k<p->gknoz; ++k)
//                     for(j=-1; j<p->gknoy; ++j)
//                         for(i=-1; i<p->gknox; ++i)
//                         {
//                             ffn=float(p->ipol4_a(eddyv));//EddyV
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                 //  Phi
//                 iin=4*(pointNum);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=-1; k<p->gknoz; ++k)
//                     for(j=-1; j<p->gknoy; ++j)
//                         for(i=-1; i<p->gknox; ++i)
//                         {
//                             ffn=float(p->ipol4phi(topo,phi));
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                 //  Elevation
//                 iin=4*(pointNum);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=0; k<p->gknoz+1; ++k)
//                     for(j=0; j<p->gknoy+1; ++j)
//                         for(i=0; i<p->gknox+1; ++i)
//                         {
//                             ffn=float(ZN[k]+(ZN[k+1]-ZN[k]));
//                             result.write((char*)&ffn, sizeof (float));
//                         }

//                 // x
//                 iin=4*(p->gknox+1);
//                 result.write((char*)&iin, sizeof (int));
//                 for(i=1; i<p->gknox+2; ++i)
//                 {
//                     ffn=float(XN[i]);
//                     result.write((char*)&ffn, sizeof (float));
//                 }
//                 // y
//                 iin=4*(p->gknoy+1);
//                 result.write((char*)&iin, sizeof (int));
//                 for(j=1; j<p->gknoy+2; ++j)
//                 {
//                     ffn=float(YN[j]);
//                     result.write((char*)&ffn, sizeof (float));
//                 }
//                 // zle
//                 iin=4*(p->gknoz+1);
//                 result.write((char*)&iin, sizeof (int));
//                 for(k=1; k<p->gknoz+2; ++k)
//                 {
//                     ffn=float(ZN[k]);
//                     result.write((char*)&ffn, sizeof (float));
//                 }

//                 result<<"\n</AppendedData>\n"
//                 <<"</VTKFile>"<<flush;

//                 result.close();
//             }
//         }
//     }
// }

void printer_CFD::print3D(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P10!=0)
    {
        if(!printerMethodInitialized)
        {
            outputMethod->setup(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed);
            printerMethodInitialized = true;
            // exit(1);
        }
        pgc->gcsync();

        std::chrono::system_clock::time_point start,end;
        if(p->mpirank==0)
            start = std::chrono::system_clock::now();
        if(outputMethod->print(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed)!=0)
            cout<<"Print out failed."<<endl;
        if(p->mpirank==0)
        {
            end = std::chrono::system_clock::now();
            auto elapsed = end - start;
            std::cout << "Print time: "<<elapsed.count() << '\n';
        }
        ++p->printcount;
        // ---------------------------------------------------------
        // comment out everything below


        // pgc->start4a(p,a->test,1);
        // pgc->start1(p,a->u,110);
        // pgc->start2(p,a->v,111);
        // pgc->start3(p,a->w,112);


        // pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
        // pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
        // pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
        // pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
        // pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
        // pgc->dgcpol4(p,a->phi,14);
        // pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
        // pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
        // pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
        // //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

        // a->u.ggcpol(p);
        // a->v.ggcpol(p);
        // a->w.ggcpol(p);
        // a->press.ggcpol(p);
        // a->eddyv.ggcpol(p);
        // a->phi.ggcpol(p);
        // a->conc.ggcpol(p);
        // a->ro.ggcpol(p);
        // a->visc.ggcpol(p);
        // a->phi.ggcpol(p);
        // a->fb.ggcpol(p);
        // a->fbh4.ggcpol(p);
        // //a->test.ggcpol(p);
        

        // pgc->gcparacox(p,a->phi,50);
        // pgc->gcparacox(p,a->phi,50);

        // pgc->gcparacox(p,a->topo,150);
        // pgc->gcparacox(p,a->topo,150);
        
        // //pgc->start4a(p,a->topo,159);

        // pgc->gcperiodicx(p,a->press,4);

        // outputFormat->extent(p,pgc);
        // if(p->mpirank==0)
        //     parallelData(a,p,pgc,pturb,pheat,pdata,pconc,pmp,psed);

        // int num=0;
        // if(p->P15==1)
        //     num = p->printcount;
        // if(p->P15==2)
        //     num = p->count;
        // int rank = p->mpirank+1;
        // outputFormat->fileName(name,"CFD",num,rank);

        // // Open File
        // ofstream result;
        // result.open(name, ios::binary);
        // if(result.is_open())
        // {
        //     n=0;

        //     offset[n]=0;
        //     ++n;

        //     // velocity
        //     offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
        //     ++n;
            
        //     pmean->offset_vtk(p,a,pgc,result,offset,n);

        //     // scalars
        //     {
        //         // pressure
        //         offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //         ++n;
        //         // k and eps
        //         pturb->offset_vtk(p,a,pgc,result,offset,n);
        //         // eddyv
        //         offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //         ++n;
        //         // phi
        //         offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //         ++n;
        //         // T
        //         pheat->offset_vtk(p,a,pgc,result,offset,n);
        //         // Multiphase
        //         pmp->offset_vtk(p,a,pgc,result,offset,n);
        //         // vorticity
        //         pvort->offset_vtk(p,a,pgc,result,offset,n);
        //         // data
        //         pdata->offset_vtk(p,a,pgc,result,offset,n);
        //         // concentration
        //         pconc->offset_vtk(p,a,pgc,result,offset,n);
        //         // rho
        //         if(p->P24==1 && p->F300==0)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // viscosity
        //         if(p->P71==1)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // VOF
        //         if(p->P72==1)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // Fi
        //         if(p->A10==4)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // conc
        //         if(p->P26==1)
        //         { 
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // topo
        //         if(p->P27==1)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // sediment bedlaod
        //         if(p->P76==1)
        //             psed->offset_vtk_bedload(p,pgc,result,offset,n);

        //         // sediment parameters 1
        //         if(p->P77==1)
        //             psed->offset_vtk_parameter1(p,pgc,result,offset,n);

        //         // sediment parameters 2
        //         if(p->P78==1)
        //             psed->offset_vtk_parameter2(p,pgc,result,offset,n);

        //         // bed shear stress
        //         if(p->P79>=1)
        //             psed->offset_vtk_bedshear(p,pgc,result,offset,n);

        //         // test
        //         if(p->P23==1)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // elevation
        //         offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //         ++n;
        //         // solid
        //         if(p->P25==1)
        //         { 
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // floating
        //         if(p->P28==1)
        //         {
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //         // walldist
        //         if(p->P29==1)
        //         {   
        //             offset[n]=offset[n-1]+4*(p->pointnum)+4;
        //             ++n;
        //         }
        //     }
        //     // end scalars
        //     outputFormat->offset(p,offset,n);
        //     //---------------------------------------------

        //     outputFormat->beginning(p,result);
            
        //     n=0;
        //     result<<"<PointData>"<<endl;
        //     result<<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //     ++n;
            
        //     pmean->name_vtk(p,a,pgc,result,offset,n);

        //     result<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //     ++n;

        //     pturb->name_vtk(p,a,pgc,result,offset,n);

        //     result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //     ++n;
        //     result<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //     ++n;

        //     pheat->name_vtk(p,a,pgc,result,offset,n);
            
        //     pmp->name_vtk(p,a,pgc,result,offset,n);

        //     pvort->name_vtk(p,a,pgc,result,offset,n);

        //     pdata->name_vtk(p,a,pgc,result,offset,n);

        //     pconc->name_vtk(p,a,pgc,result,offset,n);

        //     if(p->P24==1 && p->F300==0)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->P71==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"viscosity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }
            
        //     if(p->P72==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"VOF\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->A10==4)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->P26==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"ST_conc\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->P27==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"topo\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }
            
        //     if(p->P76==1)
        //         psed->name_vtk_bedload(p,pgc,result,offset,n);
            
        //     if(p->P77==1)
        //         psed->name_vtk_parameter1(p,pgc,result,offset,n);

        //     if(p->P78==1)
        //         psed->name_vtk_parameter2(p,pgc,result,offset,n);

        //     if(p->P79>=1)
        //         psed->name_vtk_bedshear(p,pgc,result,offset,n);

        //     if(p->P23==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     result<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //     ++n;

        //     if(p->P25==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->P28==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"floating\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }

        //     if(p->P29==1)
        //     {
        //         result<<"\t<DataArray type=\"Float32\" Name=\"walldist\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        //         ++n;
        //     }
        //     result<<"</PointData>"<<endl;
        //     outputFormat->ending(result,offset,n);
        //     //----------------------------------------------------------------------------
        //     result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

        //     //  Velocities
        //     iin=3*4*(p->pointnum);
        //     result.write((char*)&iin, sizeof (int));
        //     TPLOOP
        //     {
        //         ffn=float(p->ipol1(a->u));
        //         result.write((char*)&ffn, sizeof (float));

        //         ffn=float(p->ipol2(a->v));
        //         result.write((char*)&ffn, sizeof (float));

        //         ffn=float(p->ipol3(a->w));
        //         result.write((char*)&ffn, sizeof (float));
        //     }

        //     //  time average flow parameters
        //     pmean->print_3D(p,a,pgc,result);

        //     //  Pressure
        //     iin=4*(p->pointnum);
        //     result.write((char*)&iin, sizeof (int));
        //     TPLOOP
        //     {
        //         ffn=float(p->ipol4press(a->press)-p->pressgage);
        //         result.write((char*)&ffn, sizeof (float));
        //     }

        //     //  turbulence
        //     pturb->print_3D(p,a,pgc,result);

        //     //  eddyv
        //     iin=4*(p->pointnum);
        //     result.write((char*)&iin, sizeof (int));
        //     TPLOOP
        //     {
        //         ffn=float(p->ipol4_a(a->eddyv));
        //         result.write((char*)&ffn, sizeof (float));
        //     }

        //     //  phi
        //     nodefill4(p,a,pgc,a->phi,eta);
        //     iin=4*(p->pointnum);
        //     result.write((char*)&iin, sizeof (int));
        //     TPLOOP
        //     {
        //         if(p->P18==1)
        //             ffn=float(p->ipol4phi(a,a->phi));
        //         if(p->P18==2)
        //             ffn = float(eta(i,j,k));
        //         result.write((char*)&ffn, sizeof (float));
        //     }

        //     //  T
        //     pheat->print_3D(p,a,pgc,result);
            
        //     //  Multiphase
        //     pmp->print_3D(p,a,pgc,result);

        //     //  Vorticity
        //     pvort->print_3D(p,a,pgc,result);

        //     //  Data
        //     pdata->print_3D(p,a,pgc,result);

        //     //  Concentration
        //     pconc->print_3D(p,a,pgc,result);

        //     //  density
        //     if(p->P24==1 && p->F300==0)
        //     {
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->ro));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }
            
        //     //  viscosity
        //     if(p->P71==1)
        //     {
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4(a->visc));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }
            
        //     //  VOF
        //     if(p->P72==1)
        //     {
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4(a->vof));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }

        //     //  Fi
        //     if(p->A10==4)
        //     {
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4press(a->Fi));
        //             result.write((char*)&ffn, sizeof (float));
        //         }

        //     }

        //     if(p->P26==1)
        //     {
        //         //  conc
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4(a->conc));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }
            
        //     if(p->P27==1)
        //     {
        //         //  topo
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->topo));
        //             //ffn = float(-a->bed(i,j)+p->ZN[KP1]);
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }
            
        //     //  sediment bedload
        //     if(p->P76==1)
        //         psed->print_3D_bedload(p,pgc,result);
            
        //     //  sediment parameter 1
        //     if(p->P77==1)
        //         psed->print_3D_parameter1(p,pgc,result);

        //     //  sediment parameter 2
        //     if(p->P78==1)
        //         psed->print_3D_parameter2(p,pgc,result);

        //     //  bed shear stress
        //     if(p->P79>=1)
        //         psed->print_3D_bedshear(p,pgc,result);

        //     //  test
        //     if(p->P23==1)
        //     {
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->test));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }

        //     //  elevation
        //     iin=4*(p->pointnum)*3;
        //     result.write((char*)&iin, sizeof (int));
        //     TPLOOP
        //     {
        //         ffn=float(p->pos_z()+0.5*p->DZN[KP]);
        //         result.write((char*)&ffn, sizeof (float));
        //     }

        //     if(p->P25==1)
        //     {
        //         //  solid
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->solid));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }

        //     if(p->P28==1)
        //     {
        //         //  floating
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->fb));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }

        //     if(p->P29==1)
        //     {
        //         //  walldist
        //         iin=4*(p->pointnum);
        //         result.write((char*)&iin, sizeof (int));
        //         TPLOOP
        //         {
        //             ffn=float(p->ipol4_a(a->walld));
        //             result.write((char*)&ffn, sizeof (float));
        //         }
        //     }

        //     // -----------------------
            
        //     outputFormat->structureWrite(p,a,result);

        //     result.close();

        //     ++p->printcount;
        // }
        // else
        //     cout<<"Failed to open output file."<<endl;


        pgc->start1(p,a->u,114);
        pgc->start2(p,a->v,115);
        pgc->start3(p,a->w,116);

        // pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
        // pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
        // pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
        pgc->start4a(p,a->topo,150);

        a->u.ggcpol(p);
        a->v.ggcpol(p);
        a->w.ggcpol(p);
    }
}

// void printer_CFD::setupCompactMPIPrint(lexer *p, fdm *a, ghostcell * pgc)
// {
//     int *gorigins=nullptr;
//     int *gbeginEndPoint = nullptr;

//     double *XN = nullptr;
//     double *YN = nullptr;
//     double *ZN = nullptr;

//     int *recvcounts=nullptr;
//     int *displs=nullptr;


//     if(p->mpirank==0)
//     {
// 	    mkdir("./REEF3D_CFD_VTRCMPI",0777);

//         XN = (double *)malloc((p->gknox+2)*sizeof(double));
//         YN = (double *)malloc((p->gknoy+2)*sizeof(double));
//         ZN = (double *)malloc((p->gknoz+2)*sizeof(double));

//         gbeginEndPoint = (int *)malloc(p->mpi_size*6*sizeof(int));
//         gorigins = (int *)malloc(p->mpi_size*3*sizeof(int));

//         recvcounts = (int *)malloc(p->mpi_size*sizeof(int));
//         displs = (int *)malloc(p->mpi_size*sizeof(int));

//         cellNum=(p->gknox)*(p->gknoy)*(p->gknoz);
//         pointNum=(p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);

//         // ---------------------------------------------------------
//         // Pre-calulate offsets
//         // ---------------------------------------------------------

//         m=0;
//         compactMPIPOffset[m]=0;
//         ++m;

//         // time
//         if(p->P16==1)
//         {
//             compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(double);
//             ++m;
//         }

//         //velocities
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+3*sizeof(float)*(pointNum);
//         ++m;
//         //pressure
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
//         ++m;
//         //eddyv
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
//         ++m;
//         //phi
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
//         ++m;
//         //elevation
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
//         ++m;

//         //x
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknox+1);
//         ++m;
//         //y
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknoy+1); 
//         ++m;
//         //z
//         compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknoz+1);
//         ++m;
//     }


//     int recvcount = p->knox+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     int disp = p->origin_i+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->XN+marge,p->knox+1, XN,recvcounts,displs);

//     recvcount = p->knoy+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     disp = p->origin_j+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->YN+marge,p->knoy+1, YN,recvcounts,displs);

//     recvcount = p->knoz+1;
//     pgc->gather_int(&recvcount,1,recvcounts,1);
//     disp = p->origin_k+1;
//     pgc->gather_int(&disp,1,displs,1);
//     pgc->gatherv_double(p->ZN+marge,p->knoz+1, ZN,recvcounts,displs);

//     if(p->mpirank==0)
//     {
//         i=j=k=-1;
//         XN[0]=p->XN[IP];
//         YN[0]=p->YN[JP];
//         ZN[0]=p->ZN[KP]; 
//     }


//     kbeginPoint=-1;
//     kendPoint=p->knoz;
//     if(p->nb6>-2)
//         --kendPoint;

//     jbeginPoint=-1;
//     jendPoint=p->knoy;
//     if(p->nb2>-2)
//         --jendPoint;
    
//     ibeginPoint=-1;
//     iendPoint=p->knox;
//     if(p->nb4>-2)
//         --iendPoint;

//     int beginEndPoint[6];
//     beginEndPoint[0]=ibeginPoint;
//     beginEndPoint[1]=iendPoint;
//     beginEndPoint[2]=jbeginPoint;
//     beginEndPoint[3]=jendPoint;
//     beginEndPoint[4]=kbeginPoint;
//     beginEndPoint[5]=kendPoint;
//     pgc->gather_int(beginEndPoint,6,gbeginEndPoint,6);
//     int origins[3];
//     origins[0]=p->origin_i;
//     origins[1]=p->origin_j;
//     origins[2]=p->origin_k;
//     pgc->gather_int(origins,3,gorigins,3);

//     if(p->mpirank==0)
//     {
//         // ---------------------------------------------------------
//         // Data MPI offsets
//         // ---------------------------------------------------------
//         {
//             m=0;

//             //time
//             if(p->P16==1)
//             {
//                 offsetCMPI.push_back(compactMPIPOffset[m]);
//                 for(int n=0;n<p->mpi_size;++n)
//                 {
//                     offsetCMPIitr.push_back(offsetCMPI.size()-1);
//                     offsetCMPI.push_back(compactMPIPOffset[m]);
//                 }
//                 offsetCMPI.pop_back();
//                 m++;
//             }

//             //velocities
//             offsetCMPIPoints(p,gorigins,gbeginEndPoint,3);

//             //pressure
//             offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

//             //eddyv
//             offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

//             //phi
//             offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

//             //elevation
//             offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);
//         }

//         // header
//         {
//             header<<"<?xml version=\"1.0\"?>\n"
//             <<"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
//             <<"<RectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
//             m=0;
//             if(p->P16==1)
//             {
//                 header<<"\t<FieldData>\n";
//                 header<<"\t\t<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
//                 ++m;
//                 header<<"\t</FieldData>\n";
//             }
//             header<<"\t<Piece Extent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\">\n";
//             header<<"\t\t<PointData>\n";
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
//             ++m;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
//             ++m;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
//             ++m;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"phi\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
//             ++m;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
//             ++m;
//             header<<"\t\t</PointData>\n";
//             endIndex=m;
//             header<<"\t\t<Coordinates>\n";
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
//             m++;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
//             m++;
//             header<<"\t\t\t<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
//             m++;
//             header<<"\t\t</Coordinates>\n"
//             <<"\t</Piece>\n"
//             <<"</RectilinearGrid>\n"
//             <<"<AppendedData encoding=\"raw\">\n_";
//             headerSize=header.str().size();
//             for(auto &elem : offsetCMPI)
//                 elem+=headerSize;
//         }

//         // footer
//         {
//             // x
//             iin=4*(p->gknox+1);
//             footer.write((char*)&iin, sizeof (int));
//             for(i=1; i<p->gknox+2; ++i)
//             {
//                 ffn=float(XN[i]);
//                 footer.write((char*)&ffn, sizeof (float));
//             }
//             // y
//             iin=4*(p->gknoy+1);
//             footer.write((char*)&iin, sizeof (int));
//             for(j=1; j<p->gknoy+2; ++j)
//             {
//                 ffn=float(YN[j]);
//                 footer.write((char*)&ffn, sizeof (float));
//             }
//             // z
//             iin=4*(p->gknoz+1);
//             footer.write((char*)&iin, sizeof (int));
//             for(k=1; k<p->gknoz+2; ++k)
//             {
//                 ffn=float(ZN[k]);
//                 footer.write((char*)&ffn, sizeof (float));
//             }

//             footer<<"\n</AppendedData>\n"
//             <<"</VTKFile>";
//         }
//     }
    
//     int size=0;
//     if(p->mpirank==0)
//         size = offsetCMPI.size();
//     pgc->Bcast(&size,1,MPI_INT);
//     if(p->mpirank!=0)
//         offsetCMPI.resize(size);
//     pgc->Bcast(offsetCMPI.data(),size,MPI_OFFSET);

//     if(p->mpirank==0)
//         size = offsetCMPIitr.size();
//     pgc->Bcast(&size,1,MPI_INT);
//     if(p->mpirank!=0)
//         offsetCMPIitr.resize(size);
//     pgc->Bcast(offsetCMPIitr.data(),size,MPI_INT);
// }

// void printer_CFD::offsetCMPIPoints(lexer *p, int *gorigins, int *gbeginEndPoint, int numberOfTuples)
// {
//     offsetCMPI.push_back(compactMPIPOffset[m]);
//     for(int n=0;n<p->mpi_size;++n)
//     {
//         offsetCMPIitr.push_back(offsetCMPI.size()-1);
//         if(n>0)
//         {
//             offsetCMPI.pop_back();
//             offsetCMPI.push_back(compactMPIPOffset[m]+(gorigins[0+3*n]+gorigins[1+3*n]*(p->gknox+1)+gorigins[2+3*n]*(p->gknox+1)*(p->gknoy+1))*numberOfTuples*sizeof(float)+sizeof(int));
//         }
//         for(int k=0;k<(gbeginEndPoint[5+6*n]-gbeginEndPoint[4+6*n]);k++)
//         {
//             if(k>0)
//             {
//                 offsetCMPI.pop_back();
//                 offsetCMPI.push_back(compactMPIPOffset[m]+(gorigins[0+3*n]+gorigins[1+3*n]*(p->gknox+1)+(gorigins[2+3*n]+k)*(p->gknox+1)*(p->gknoy+1))*numberOfTuples*sizeof(float)+sizeof(int));
//             }
//             for(int j=0;j<(gbeginEndPoint[3+6*n]-gbeginEndPoint[2+6*n]);j++)
//             {
//                 offsetCMPI.push_back(offsetCMPI.back()+numberOfTuples*sizeof(float)*(p->gknox+1));
//                 if(n==0&&k==0&&j==0)
//                     offsetCMPI.back()+=sizeof(int);
//             }
//         }
//     }
//     offsetCMPI.pop_back();
//     m++;
// }

// void printer_CFD::print3DcompactMPI(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
// {
//     if(p->P10==2)
//     {
//         pgc->start4a(p,a->test,1);
//         pgc->start1(p,a->u,110);
//         pgc->start2(p,a->v,111);
//         pgc->start3(p,a->w,112);


//         pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
//         pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
//         pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
//         pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
//         pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
//         pgc->dgcpol4(p,a->phi,14);
//         pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
//         pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
//         pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
//         //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

//         a->u.ggcpol(p);
//         a->v.ggcpol(p);
//         a->w.ggcpol(p);
//         a->press.ggcpol(p);
//         a->eddyv.ggcpol(p);
//         a->phi.ggcpol(p);
//         a->conc.ggcpol(p);
//         a->ro.ggcpol(p);
//         a->visc.ggcpol(p);
//         a->phi.ggcpol(p);
//         a->fb.ggcpol(p);
//         a->fbh4.ggcpol(p);
//         //a->test.ggcpol(p);
        

//         pgc->gcparacox(p,a->phi,50);
//         pgc->gcparacox(p,a->phi,50);

//         pgc->gcparacox(p,a->topo,150);
//         pgc->gcparacox(p,a->topo,150);
        
//         //pgc->start4a(p,a->topo,159);

//         pgc->gcperiodicx(p,a->press,4);


//         MPI_File file;
//         if(p->mpirank==0)
//         {
//             int num=0;
//             if(p->P15==1)
//                 num = p->printcount-1;
//             if(p->P15==2)
//                 num = p->count;
//             snprintf(name,200,"./REEF3D_CFD_VTRCMPI/REEF3D-CFD-%08i.vtr",num);
//         }
//         pgc->Bcast(name,200,MPI_CHAR);
//         pgc->File_open_createWriteOnly(&file,name);
//         pgc->File_set_size(file,0);

//         // header
//         if(p->mpirank==0)
//         {
//             pgc->File_write_at_char(file, 0, header.str().c_str(), header.str().size());
//         }

//         // data + footer
//         {
//             std::stringstream result;
//             m=0;

//             // Time
//             if(p->P16==1)
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=sizeof(double);
//                     result.write((char*)&iin, sizeof (int));
//                     std::stringstream time;
//                     time<<std::setprecision(7)<<p->simtime;
//                     double t = std::stod(time.str());
//                     result.write((char*)&t, sizeof(double));
//                     pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m]], result.str().c_str(), result.str().size());
//                     result.str(std::string());
//                     result.clear();
//                 }
//                 ++m;
//             }

//             //  Velocities
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=3*sizeof(float)*(pointNum);
//                     result.write((char*)&iin, sizeof (int));
//                 }
//                 int n=0;
//                 for(k=kbeginPoint;k<kendPoint;++k)
//                     for(j=jbeginPoint;j<jendPoint;++j)
//                     {
//                         for(i=ibeginPoint;i<iendPoint;++i)
//                         {
//                             ffn=float(p->ipol1(a->u));
//                             result.write((char*)&ffn, sizeof (float));

//                             ffn=float(p->ipol2(a->v));
//                             result.write((char*)&ffn, sizeof (float));

//                             ffn=float(p->ipol3(a->w));
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                         pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
//                         result.str(std::string());
//                         result.clear();
//                         ++n;
//                     }
//                 ++m;
//             }

//             //  Pressure
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=3*sizeof(float)*(pointNum);
//                     result.write((char*)&iin, sizeof (int));
//                 }
//                 int n=0;
//                 for(k=kbeginPoint;k<kendPoint;++k)
//                     for(j=jbeginPoint;j<jendPoint;++j)
//                     {
//                         for(i=ibeginPoint;i<iendPoint;++i)
//                         {
//                             ffn=float(p->ipol4press(a->press)-p->pressgage);
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                         pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
//                         result.str(std::string());
//                         result.clear();
//                         ++n;
//                     }
//                 ++m;
//             }

//             //  EddyV
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=sizeof(float)*(pointNum);
//                     result.write((char*)&iin, sizeof (int));
//                 }
//                 int n=0;
//                 for(k=kbeginPoint;k<kendPoint;++k)
//                     for(j=jbeginPoint;j<jendPoint;++j)
//                     {
//                         for(i=ibeginPoint;i<iendPoint;++i)
//                         {
//                             ffn=float(p->ipol4_a(a->eddyv));
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                         pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
//                         result.str(std::string());
//                         result.clear();
//                         ++n;
//                     }
//                 ++m;
//             }

//             //  Phi
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=sizeof(float)*(pointNum);
//                     result.write((char*)&iin, sizeof (int));
//                 }
//                 int n=0;
//                 nodefill4(p,a,pgc,a->phi,eta);
//                 for(k=kbeginPoint;k<kendPoint;++k)
//                     for(j=jbeginPoint;j<jendPoint;++j)
//                     {
//                         for(i=ibeginPoint;i<iendPoint;++i)
//                         {
//                             if(p->P18==1)
//                                 ffn=float(p->ipol4phi(a,a->phi));
//                             if(p->P18==2)
//                                 ffn = float(eta(i,j,k));
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                         pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
//                         result.str(std::string());
//                         result.clear();
//                         ++n;
//                     }
//                 ++m;
//             }

//             //  Elevation
//             {
//                 if(p->mpirank==0)
//                 {
//                     iin=sizeof(float)*(pointNum);
//                     result.write((char*)&iin, sizeof (int));
//                 }
//                 int n=0;
//                 for(k=kbeginPoint;k<kendPoint;++k)
//                     for(j=jbeginPoint;j<jendPoint;++j)
//                     {
//                         for(i=ibeginPoint;i<iendPoint;++i)
//                         {
//                             ffn=float(p->pos_z()+0.5*p->DZN[KP]);
//                             result.write((char*)&ffn, sizeof (float));
//                         }
//                         pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
//                         result.str(std::string());
//                         result.clear();
//                         ++n;
//                     }
//                 ++m;
//             }

//             // footer
//             if(p->mpirank==0)
//             {
//                 pgc->File_write_at_char(file, headerSize + compactMPIPOffset[m], footer.str().c_str(), footer.str().size());
//             }
//         }
//         pgc->File_close(&file);
//     }
// }
