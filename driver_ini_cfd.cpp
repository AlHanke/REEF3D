/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"waves_header.h"
#include"lexer.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"
#include<sys/stat.h>
#include<sys/types.h>

void driver::driver_ini()
{
p->count=0;

p->cellnumtot=pgc->globalisum(p->cellnum);
p->pointnumtot=pgc->globalisum(p->pointnum);

if(p->mpirank==0)
cout<<"number of cells: "<<p->cellnumtot<<endl;

	log_ini();

if(p->mpirank==0)
cout<<"starting driver_ini"<<endl;
    
    
	// Solid
    if(p->G39==1)
    {
    solid solid_object(p,a,pgc);
    solid_object.start(p,a,pgc,pflow,pconvec,preso);
    }
    
    // Geotopo
    if((p->G50>0 && p->G51>0) || p->G60>0 || p->G61>0)
    {
    geotopo gtopo(p,a,pgc);
    gtopo.start(p,a,pgc,pflow,pconvec,preto,pvrans);
    }
    
	// 6DOF
	if(p->X10==1 && p->X13!=2)
    {}
    else
    {
	    p6dof->initialize(p,a,pgc,pnet);
    }

    // Sediment
	if(p->S10>0)
    {
    psed->ini(p,a,pgc);
    for(int qn=0;qn<5;++qn)
    psed->relax(p,a,pgc);
    preto->start(a,p,a->topo,pconvec,pgc);
    psed->update(p,a,pgc,pflow);
    pgc->start4a(p,a->topo,150);
    }
    
    //ioflow ini
    pflow->ini(p,a,pgc);

    
	starttime=pgc->timer();
    if(p->B60>0 || p->T36==2)
	pgc->walldistance(p,a,pgc,pconvec,preini,pflow,a->walld);
	
	pflow->inflow_walldist(p,a,pgc,pconvec,preini,pflow);
	
	double walltime=pgc->timer()-starttime;
	
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"Walldist time: "<<setprecision(4)<<walltime<<endl;

	if(p->P150==0)
	pdata = new data_void(p,a,pgc);
	
	if(p->P150>0)
	pdata = new data_f(p,a,pgc);
	
	pdata->start(p,a,pgc);

    pnse->ini(p,a,pgc,pflow);
	
    pheat->heat_ini(p,a,pgc,pheat);
	pconc->ini(p,a,pgc,pconc);

    ptstep->ini(a,p,pgc);
    pini->iniphi_io(a,p,pgc);
	pflow->gcio_update(p,a,pgc);
	pflow->pressure_io(p,a,pgc);
    
    if (p->F80>0)
    {
        pflow->vof_relax(p,pgc,a->vof);
    }
    
    else
    if(p->F30>0 || p->F40>0)
    {
        for(int qn=0;qn<20;++qn)
        pflow->phi_relax(p,pgc,a->phi);
        preini->start(a,p, a->phi, pgc, pflow);
        pfsf->update(p,a,pgc,a->phi);        
        pini->iniphi_surfarea(p,a,pgc);
    }
    
	ppart->setup(p,a,pgc);
	pini->iniphi_io(a,p,pgc);
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	potflow->start(p,a,ppoissonsolv,pgc);
    pflow->wavegen_precalc(p,pgc);
	if(p->I12>=1)
	pini->hydrostatic(p,a,pgc);

	if(p->I11==1)
	ptstep->start(a,p,pgc,pturb);
    
    if(p->I13==1)
    pturb->ini(p,a,pgc);
	

    if(p->I58_2>0.0)
	pini->droplet_ini(p,a,pgc);

	pflow->pressure_io(p,a,pgc);
    
    poneph->update(p,a,pgc,pflow);

	pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);

    pgc->start4(p,a->press,40);

	if(p->I40==1)
	pini->stateini(p,a,pgc,pturb);
    
	pgc->start4(p,a->press,40);
	

    pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,psed);

// ini variables
    for(int qn=0; qn<2; ++qn)
    {
    pturb->ktimesave(p,a,pgc);
    pturb->etimesave(p,a,pgc);
    }

    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

}

