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
#include <string>

/// @brief Sediment model on particle basis
/// Class handling the sediment when using the options for Lagrangian particles and VRANS.\n
/// Does all the initialization of the topography with particles, modification of topo values and print out.
/// @param p 
/// @param pgc 
/// @param pturb 
sedpart::sedpart(lexer* p, ghostcell* pgc, turbulence *pturb) : particle_func(p), PP(10,p->S20,p->S22,true), active_box(p), active_topo(p), irand(10000), drand(irand)
{
    pvrans = new vrans_f(p,pgc);
    pbedshear  = new bedshear(p,pturb);
    printcount = 0;

    // Create Folder
	if(p->mpirank==0 && (p->Q180>0||p->Q182>0))
	    mkdir("./REEF3D_CFD_SedPart",0777);

    // Output configuration to console
    if(p->mpirank==0)
    {
        string buff;
        buff.append("\nSedPart configuration\nParticles and VRANS active\n");
        buff.append("General configuration:\n\tTopo deformation: ");
        p->Q13>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tBox seeding: ");
        p->Q110>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tPoint seeding: ");
        p->Q61>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tTopo seeding: ");
        p->Q101>0?buff.append("True\n"):buff.append("False\n");
        buff.append("\tSuspension seeding: ");
        p->Q120>0?buff.append("True\n"):buff.append("False\n");
        buff.append("Particle properties:\n\td50: "+std::to_string(p->S20)+" m\n\tDensity: "+std::to_string(p->S22)+" kg/m/m/m\n\tPorosity: "+std::to_string(p->S24)+"\n");
        buff.append("Seeding properties:\n\tSeed: "+(p->Q29>0?std::to_string(p->Q29):"time dep.")+"\n\tParticles per cell: "+std::to_string(p->Q24)+"\n\tParticles represened by one: "+std::to_string(p->Q41)+"\n");
        cout<<buff<<endl;
    }
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
        /// runtime seeding
		if(p->Q120==1&&p->count%p->Q121==0)
			posseed_suspended(p,a);
        point_source(p,a);

        /// transport
        // pgc->gcparaxijk(p,cellSum,1); ghostcell exchange needed
        erode(p,a,pgc);
        transport(p,a,&PP);
		xchange=transfer(p,pgc,&PP,maxparticle);
		removed=remove(p,&PP);

        /// topo update
		make_stationary(p,a,&PP);
        particlesPerCell(p,pgc,&PP);
        if(p->Q13==1)
        update_cfd(p,a,pgc,pflow,preto);
        particleStressTensor(p,a,pgc,&PP);
        /// cleanup
        if(p->count%p->Q20==0)
        {
            if(PP.size == 0)
                PP.erase_all();
            // PP.optimize();
            cleanup(p,a,&PP,0);
        }
	}

    /// print out
	print_particles(p,a,pgc);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);
	p->sedsimtime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<"\nTotal bed volume change: "<<std::setprecision(9)<<volumeChangeTotal<<endl;

    // testing
}

void sedpart::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    // seed
    seed_ini(p,a,pgc);
    gpartnum=pgc->globalisum(partnum);
    allocate(p);
    seed(p,a);
    make_stationary(p,a,&PP);
    particlesPerCell(p,pgc,&PP);
    particleStressTensor(p,a,pgc,&PP);
    
    // print
    print_vtu(p,a,pgc);
    printcount++;
    gparticle_active = pgc->globalisum(PP.size);
    if(p->mpirank==0)
        cout<<"Sediment particles: active: "<<gparticle_active<<endl;
    
    // vrans
    pvrans->sed_update(p,a,pgc);

    // testing
    PLAINLOOP
    a->test(i,j,k)=active_topo(i,j,k);
    volumeChangeTotal=0;
    if(0==p->mpirank)
    {
    }
    // ILOOP
    // if(p->XN[IP]==0.2525)
    // cout<<p->mpirank<<endl;

    if(p->mpirank==2)
    {
        // ILOOP
        // for(int q=0;q<p->margin;++q)
        // cout<<"Topo after ini("<<i<<"): "<<a->topo(i,-q,19)<<"|"<<a->topo(i,-q,20)<<endl;
        // int qq;
        // QQGC4A
        // if(p->gcb4a[qq][0]==32&&p->gcb4a[qq][1]==0&&p->gcb4a[qq][2]==20)
        for(int n=0;n<6;n++)
        cout<<a->topo(31+n,0,20)<<endl;
    }
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
    ILOOP
    JLOOP
    if(p->flag_topo_changed[IJ]==1)
    {
        double dh=p->topo_change[IJ]/p->DXN[IP]/p->DYN[JP];
        a->bed(i,j)+=dh;
        KLOOP
        {
            a->topo(i,j,k) -= dh;
            // Seeding update
            if((abs(a->topo(i,j,k))<(p->DZN[KP]*ceil(p->Q102)))&&(a->topo(i,j,k)<=0.25*p->DZN[KP]))
                active_topo(i,j,k) = 1.0;
            else
                active_topo(i,j,k) = 0.0;
        }
        volumeChangeTotal += p->topo_change[IJ];

        // Reset
        p->flag_topo_changed[IJ]=0;
        p->topo_change[IJ]=0;
    }

    pgc->start4a(p,a->topo,150);
    pgc->gcsl_start4(p,a->bed,50);
    preto->start(p,a,pgc,a->topo);
    if(p->mpirank==0)
        cout<<"Topo: update grid..."<<endl;
    pvrans->sed_update(p,a,pgc);
    pflow->gcio_update(p,a,pgc);

    // fixPos
    // fixPos(p,a,&PP);
    // if(p->Q101>0)
    // posseed_topo(p,a);
}

void sedpart::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
}

void sedpart::erode(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->Q101>0)
        make_moving(p,a,&PP);
    
    
    // int i,j,k;
    // double x,y,z;
    // double eroded=0.0;
    // size_t index=0;
    // // pbedshear->taueff_loc
    // cout<<"eroding..."<<endl;
    // if(p->count%p->Q121==0)
    //     SLICEBASELOOP
    //     {
    //         // test for erosion
    //         if (i%2==0&&j%3==0)
    //             eroded +=volume(&PP,0);

    //         // Change amount accoding to eroded volume?
    //         // Rerun eroded column until no more erosion?

    //         while (eroded>0)
    //         {
    //             x = p->XN[IP] + p->DXN[IP]*double(rand() % irand)/drand;
    //             y = p->YN[JP] + p->DYN[JP]*double(rand() % irand)/drand;
    //             z = p->ZN[KP] + p->DZN[KP]*0.5;
    //             z-=p->ccipol4_b(a->topo,x,y,z);
    //             cout<<PP.size<<"|";
    //             index=PP.add(x,y,z,1);
    //             cout<<PP.size<<endl;
    //             eroded -= volume(&PP,index);
    //         }
    //         eroded=0;
    //     }
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