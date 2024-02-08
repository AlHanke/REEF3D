/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "sedpart.h"

#include "lexer.h"
#include "ghostcell.h"
#include "looping.h"
#include "fdm.h"
#include "reinitopo.h"
#include "vrans_f.h"
#include "ioflow.h"
#include "turbulence.h"
#include "bedshear.h"

#include <sys/stat.h>

sedpart::sedpart(lexer* p, ghostcell* pgc, turbulence *pturb) : particle_func(p), PP(10,p->S20,p->S22,p->S24), active_box(p), active_topo(p), irand(10000), drand(irand)
{
    pvrans = new vrans_f(p,pgc);
    pbedshear  = new bedshear(p,pturb);
    PP.ini_cellSum(p->imax*p->jmax*p->kmax);

    // Create Folder
	if(p->mpirank==0 && p->P14==1 && p->Q180>0)
	    mkdir("./REEF3D_CFD_SedPart",0777);
}
sedpart::~sedpart()
{

}

/// @brief 
///
/// @param p 
/// @param a 
/// @param pgc 
/// @param pflow 
/// @param preto 
/// @param psolv 
void sedpart::start_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow,
                                    reinitopo* preto, solver* psolv)
{
    double starttime=pgc->timer();
	int xchange=0;
	int removed=0;

	if (p->count>=p->Q43)
	{
		if(p->Q120==1&&p->count%p->Q121==0)
			posseed_suspended(p,a,pgc);
        // erode(p,a,pgc);
		// advect(p,a,&PP,0,0,0,0);
        transport(p,a,&PP);
		xchange=transfer(p,pgc,&PP,maxparticle);
		removed=remove(p,&PP);
		make_stationary(p,a,&PP);
        if(p->Q13==1)
            update_cfd(p,a,pgc,pflow,preto);
	}

	print_particles(p,a,pgc);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);
	p->sedsimtime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: active: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sed. part. sim. time: "<<p->sedsimtime<<endl;
}

void sedpart::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    // seed
    seed_ini(p,a,pgc);
    gpartnum=pgc->globalisum(partnum);
    allocate(p,a,pgc);
    seed(p,a,pgc);
    make_stationary(p,a,&PP);
    
    // print
    print_vtu(p,a,pgc);
    printcount++;
    gparticle_active = pgc->globalisum(PP.size);
    if(p->mpirank==0)
        cout<<"Sediment particles: active: "<<gparticle_active<<endl;
}

void sedpart::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow*, slice &P, slice &Q)
{
}

void sedpart::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
}
    
void sedpart::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo* preto)
{
    int i,j,k;
    JILOOP
        if(p->flag_topo_changed[IJ]==1)
        {
            double dh=p->topo_change[IJ]/p->DXN[IP]/p->DYN[JP]*100;
            cout<<"dh["<<p->mpirank<<"]("<<i<<","<<j<<")="<<dh<<"|"<<IJ<<endl;
            KLOOP
                a->topo(i,j,k) -= dh;
            p->flag_topo_changed[IJ]=0;
            p->topo_change[IJ]=0;
        }

    pgc->start4a(p,a->topo,150);
    preto->start(p,a,pgc,a->topo);
    if(p->mpirank==0)
        cout<<"Topo: update grid..."<<endl;
    pvrans->sed_update(p,a,pgc);
    pflow->gcio_update(p,a,pgc);
}

void sedpart::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
}

void sedpart::erode(lexer* p, fdm* a, ghostcell* pgc)
{
    int i,j,k;
    double x,y,z;
    double eroded=0.0;
    size_t index=0;
    // pbedshear->taueff_loc
    cout<<"eroding..."<<endl;
    if(p->count%p->Q121==0)
        SLICEBASELOOP
        {
            // test for erosion
            if (i%2==0&&j%3==0)
                eroded +=volume(&PP,0);

            // Change amount accoding to eroded volume?
            // Rerun eroded column until no more erosion?

            while (eroded>0)
            {
                x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
                y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
                z = p->ZN[KP] + p->DZN[KP]*0.5;
                z-=p->ccipol4_b(a->topo,x,y,z);
                cout<<PP.size<<"|";
                index=PP.add(x,y,z,1);
                cout<<PP.size<<endl;
                eroded -= volume(&PP,index);
            }
            eroded=0;
        }
}

void sedpart::relax(lexer *p,ghostcell *pgc)
{
}

double sedpart::qbeval(int ii, int jj)
{
    double val=0.0;

    return val;
}

void sedpart::qbeget(int ii, int jj, double val)
{
}

double sedpart::bedzhval(int ii, int jj)
{
    double val=0.0;

    return val;
}

void sedpart::print_2D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::print_3D_bedload(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::name_pvtu_bedload(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sedpart::name_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtp_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtu_bedload(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::print_2D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::print_3D_bedshear(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::name_pvtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sedpart::name_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtp_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtu_bedshear(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::print_2D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::print_3D_parameter1(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::name_pvtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sedpart::name_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtp_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtu_parameter1(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::print_2D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::print_3D_parameter2(lexer* p, ghostcell *pgc, ofstream &result)
{	
}

void sedpart::name_pvtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result)
{
}

void sedpart::name_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtp_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sedpart::offset_vtu_parameter2(lexer *p, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

double sedpart::bedshear_point(lexer *p, fdm *a,ghostcell *pgc)
{
	return 0.0;
}

void sedpart::start_susp(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
}

void sedpart::ctimesave(lexer *p, fdm* a)
{

}