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
#include "sedpart_movement.h"

#include "lexer.h"
#include "ghostcell.h"
#include "fdm.h"
#include "vrans_f.h"
#include "reinitopo.h"
#include "ioflow.h"
#include "bedshear.h"
#include "sediment_fdm.h"

/// @brief Initializes everything in the sediment for the CFD solver
/// Determines cell which should be filled with particles
/// Allocates memory for the particles
/// Seeds the particles
/// Prepares particles and particle related variables for simulation
/// Initializes VRANS
void sedpart::ini_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    // vrans
    pvrans->sed_update(p,a,pgc);
    ALOOP
	{
        if(a->topo(i,j,k)<0.0)
	        a->porosity(i,j,k)= p->S24;
        else
	        a->porosity(i,j,k)=1.0;
    }

    movement->setup(p,*a,PP.d50);

    volume0 = movement->volume(p,*a,PP);
    volume0 = pgc->globalsum(volume0);
    volume = volume0;

    if(p->I40!=1)
    {
        // seed
        seed_ini(p,a,pgc);
        PP.reserve(maxparticle);
        seed(p,a);
        make_stationary(p,a,&PP);
    }

    gparticle_active = pgc->globalisum(PP.size);

    fill_PQ_cfd(p,a,pgc);
    // movement->move(p,*a,*pgc,PP,s,*pturb);
    if(p->Q13==1)
        movement->update(p,*a,*pgc,PP);
    
    // print
    if((p->I40!=1)||(p->I40==1&&inicount>0))
    print_particles(p);
    
    if(p->mpirank==0)
    {
        if(p->I40!=1)
            cout<<"Sediment particles: "<<gparticle_active<<"\n";
        else if (inicount>0)
            cout<<"Loaded particles "<<gparticle_active<<" from state file.\n";
        cout<<"Initial bed volume: "<<std::setprecision(prec)<<volume0<<" m^3"<<endl;
    }
    
    SLICELOOP4
    s.bedk(i,j)=0;
    
    SLICELOOP4
    {
        KLOOP
            PBASECHECK
            if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
                s.bedk(i,j)=k+1;
        s.reduce(i,j)=0.3;
    }
    pbedshear.taubed(p,a,pgc,&s);
    pgc->gcsl_start4(p,s.tau_eff,1);
    pbedshear.taucritbed(p,a,pgc,&s);
    pgc->gcsl_start4(p,s.tau_crit,1);
    
    ++inicount;
    debug(p,a,pgc);

    // if(inicount==1)
    // {
    //     particles_obj Dummy(1000);
    //     seedDummy(p,a,Dummy);
    //     printDummyVTP(p,Dummy);
    // }
}

/// @brief CFD calculation function
/// @param a fdm object
/// @param pflow IO-flow object
/// @param preto topography reinitialization object
void sedpart::start_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow,
                                    reinitopo* preto, solver* psolv)
{
    double starttime=pgc->timer();
	int xchange=0;
	int removed=0;

    if (p->count>=p->Q43)
	{
        fill_PQ_cfd(p,a,pgc);
        pslope.slope_cds(p,pgc,&s);
        pbedshear.taubed(p,a,pgc,&s);
        preduce->start(p,pgc,&s);
        pgc->gcsl_start4(p,s.tau_eff,1);
        pbedshear.taucritbed(p,a,pgc,&s);
        pgc->gcsl_start4(p,s.tau_crit,1);
	
        /// runtime seeding
		if(p->Q120==1&&p->count%p->Q121==0)
			posseed_suspended(p,a);
        point_source(p,a);
        if(p->Q101>0)
        {
            // topo_influx(p,a);
            // solid_influx(p,a);
            // set_active_topo(p,a);
            // posseed_topo(p,a);
            // seed_topo(p,a);
        }

        /// transport
        erode(p,a);
        movement->move(p,*a,*pgc,PP,s,*pturb);
		xchange=transfer(p,pgc,&PP, *movement, maxparticle);
		removed=remove(p,&PP);
        deposit(p,a);

        /// topo update
        update_cfd(p,a,pgc,pflow,preto);

        /// cleanup
        if(p->Q20>=0 && p->count%p->Q20==0)
        {
            if(PP.size == 0)
                PP.erase_all();
            PP.optimize();
        }
	}

    /// print out
	print_particles(p);

	gparticle_active = pgc->globalisum(PP.size);
    gremoved = pgc->globalisum(removed);
    gxchange = pgc->globalisum(xchange);

    volume = movement->volume(p,*a,PP);
    volume = pgc->globalsum(volume);

	p->sedsimtime=pgc->timer()-starttime;

    if(p->mpirank==0 && (p->count%p->P12==0))
    	cout<<"Sediment particles: "<<gparticle_active<<" | xch: "<<gxchange<<" rem: "<<gremoved<<" | sim. time: "<<p->sedsimtime<<"\nTotal bed volume change: "<<(volume-volume0)/volume0<<" %, "<<std::setprecision(prec)<<volume-volume0<<" m^3"<<endl;
    debug(p,a,pgc);
}

/// @brief Updates the topography for the CFD solver
void sedpart::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo* preto)
{
    prelax.start(p,pgc,&s);
    if(p->Q13==1)
        movement->update(p,*a,*pgc,PP);
    preto->start(p,a,pgc,a->topo);
    if(p->mpirank==0)
        cout<<"Topo: update grid..."<<endl;
    prelax.start(p,pgc,&s);
    pvrans->sed_update(p,a,pgc);
    pflow->gcio_update(p,a,pgc);
}

void sedpart::fill_PQ_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    double zval,xip,yip;

    SLICELOOP4
        s.bedk(i,j)=0;
    
    SLICELOOP4
        KLOOP
            PBASECHECK
                if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
                    s.bedk(i,j)=k+1;
    
    SLICELOOP1
    {
        k=s.bedk(i,j);
        
        xip= p->XN[IP1];
        yip= p->YP[JP];
        zval = 0.5*(s.bedzh(i,j)+s.bedzh(i+1,j)) + 1.6*p->DZN[k];
        
        s.P(i,j) = a->P(i,j) = p->ccipol1_a(a->u,xip,yip,zval);
    }
    
    SLICELOOP2
    {
        k=s.bedk(i,j);
        
        xip= p->XP[IP];
        yip= p->YN[JP1];
        zval = 0.5*(s.bedzh(i,j)+s.bedzh(i,j+1)) + 1.6*p->DZN[k];
        
        s.Q(i,j) = a->Q(i,j)  = p->ccipol2_a(a->v,xip,yip,zval);
    }
    
    pgc->gcsl_start1(p,s.P,10);
	pgc->gcsl_start2(p,s.Q,11);
    
    pgc->gcsl_start1(p,a->P,10);
	pgc->gcsl_start2(p,a->Q,11);
}