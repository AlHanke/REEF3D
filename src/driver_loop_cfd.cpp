/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include <cstdlib>
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"6DOF_header.h"

void driver::loop_cfd(fdm* a)
{
    if(p->mpirank==0)
    cout<<"starting mainloop.CFD"<<endl;

    // Per-stage velocity growth probe (REEF_STEP_PROBE): max|u/v/w| over valid cells
    // per level at named points in the step, to localise which stage first grows the
    // velocity after the init regrid (predictor vs projection vs free surface).
    auto stepprobe = [&](const char* tag)
    {
        #if USE_AMREX
        if(!std::getenv("REEF_STEP_PROBE")) return;
        double mu=0,mv=0,mw=0;
        for(int lev=0; lev<p->nlevs; ++lev)
        {
            mu = std::max(mu, a->u.GetMultiFab(lev).norm0());
            mv = std::max(mv, a->v.GetMultiFab(lev).norm0());
            mw = std::max(mw, a->w.GetMultiFab(lev).norm0());
        }
        if(p->mpirank==0)
            std::cout<<"  [stepprobe] count="<<p->count<<" "<<tag
                     <<"  |u|="<<mu<<" |v|="<<mv<<" |w|="<<mw<<std::endl;
        #endif
    };
    
    //vec_test(p,a,pgc,a->test);
    //pos_test(p,a,pgc);
    //ipol_test(p,a,pgc);
    //bedslope_test(p,pgc);
    
//-----------MAINLOOP CFD----------------------------
    while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
    {
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
            cout<<"------------------------------------"<<endl;
            cout<<p->count<<endl;

            cout<<"simtime: "<<p->simtime<<endl;
            cout<<setprecision(5)<<"timestep: "<<p->dt<<endl;

            if(p->X10>0)
            cout<<"fbtimestep: "<<p->fbdt<<" fbmax: "<<p->fbmax<<endl;

            if(p->B90>0 && p->B92<=11)
            cout<<"t/T: "<<p->simtime/p->wT<<endl;

            if(p->B90>0 && p->B92>11)
            cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }
        double temp_time0a = pgc->timer();
        pflow->flowfile(p,a,pgc,pturb);

        pflow->wavegen_precalc(p,pgc);


        double temp_time0c = pgc->timer();
        pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,a->phi);
        double temp_time0d = pgc->timer();
        pturb->start(a,p,pturbdisc,pturbdiff,psolv,pgc,pflow,pvrans);
        pheat->start(a,p,pheatdisc,pdiff,psolv,pgc,pflow);
        pconc->start(a,p,pconcdisc,pconcdiff,pturb,psolv,pgc,pflow);
        pmp->start(p,a,pgc,pmpconvec,psolv,pflow,preini);

        psed->start_susp(p,a,pgc,pflow,psolv);
        psed->start_cfd(p,a,pgc,pflow,preto,psolv);
        double temp_time0b = pgc->timer();
        
        double temp_time1 = pgc->timer();
        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        double flow_relax_time = pgc->timer();
        stepprobe("pre-fsf");
        pfsf->update(p,a,pgc,a->phi);
        stepprobe("post-fsf");
        double fsf_update_time = pgc->timer();
        // B5: re-equilibrate press against the post-reinit density before the predictor. Must use
        // ppoissonsolv (the pressure/hypre_ssamg solver that owns cf_links/cf_masks/matvec), NOT
        // the generic psolv -- the projection's C-F coupling and solve group 5 live there.
        ppress->rebalance(p,a,pgc,ppois,ppoissonsolv,pflow);

        // EXPERIMENT (REEF_HYDRO_REBUILD): the level set is re-distanced in pfsf->update above, so
        // roface changes while press still holds last step's value -> the predictor hydrostatic
        // gradient is stale -> a per-step surface seed the multi-level coupling amplifies. Rebuild
        // the hydrostatic press against the current phi so grad(press)=W22*roface here. NOTE: wipes
        // dynamic pressure -> valid only for the at-rest hydrostatic test; the production fix must
        // split press into a rebuilt hydrostatic reference + solved dynamic part.
        if(std::getenv("REEF_HYDRO_REBUILD")) pini->hydrostatic(p,a,pgc);
        stepprobe("pre-mom");
        pmom->start(p,a,pgc,pvrans,p6dof);
        stepprobe("post-mom");
        double momentum_time = pgc->timer();
        pbench->start(p,a,pgc,pconvec);

        //save previous timestep
        double temp_time2 = pgc->timer();
        pturb->ktimesave(p);
        pturb->etimesave(p);
        pflow->veltimesave(p,a,pgc,pvrans);
        psed->ctimesave(p,a);
        double save_time = pgc->timer();

        //time advancement
        p->simtime+=p->dt;

        // printer
        double temp_time3 = pgc->timer();
        pprint->start(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);
        double print_time = pgc->timer();

        p->regrid(a,preini,p6dof,pgc,pflow);
        // ppress->rebalance(p,a,pgc,ppois,psolv,pflow);

        //timestep control
        ptstep->start(a,p,pgc,pturb);

        // Shell-Printout
        if(p->mpirank==0)
        {
            endtime=pgc->timer();

            p->itertime=endtime-starttime;
            p->totaltime+=p->itertime;
            p->gctotaltime+=p->gctime;
            p->Xtotaltime+=p->xtime;
            p->meantime=(p->totaltime/double(p->count));
            p->gcmeantime=(p->gctotaltime/double(p->count));
            p->Xmeantime=(p->Xtotaltime/double(p->count));

            if(p->count%p->P12==0)
            {
                if(p->B90>0)
                cout<<"wavegentime: "<<setprecision(3)<<p->wavecalctime<<endl;
                if(p->X10>0)
                cout<<"fbtime: "<<setprecision(3)<<p->fbtime<<endl;
                cout<<"reinitime: "<<setprecision(3)<<p->reinitime<<endl;
                cout<<"gctime: "<<setprecision(3)<<p->gctime<<"\t average gctime: "<<setprecision(3)<<p->gcmeantime<<endl;
                cout<<"Xtime: "<<setprecision(3)<<p->xtime<<"\t average Xtime: "<<setprecision(3)<<p->Xmeantime<<endl;
                cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
                cout<<"timer per step: "<<setprecision(3)<<p->itertime<<endl;

                if(std::getenv("REEF_timing"))
                {
                    const int precision = 6;
                    const double t1 = temp_time0b-temp_time0a;
                    const double t2 = flow_relax_time-temp_time1;
                    const double t3 = fsf_update_time-flow_relax_time;
                    const double t4 = momentum_time-fsf_update_time;
                    const double t5 = save_time-temp_time2;
                    const double t6 = print_time-temp_time3;
                    const double t7 = temp_time0d-temp_time0c;
                    const double iter_time_woutprint = p->itertime - t6;
                    cout<<"precalc time: "<<setprecision(precision)<<t1<<":"<<t1/iter_time_woutprint*100<<endl;
                    cout<<"fsf start time: "<<setprecision(precision)<<t7<<":"<<t7/iter_time_woutprint*100<<endl;
                    cout<<"flow relax time: "<<setprecision(precision)<<t2<<":"<<t2/iter_time_woutprint*100<<endl;
                    cout<<"fsf update time: "<<setprecision(precision)<<t3<<":"<<t3/iter_time_woutprint*100<<endl;
                    cout<<"momentum time: "<<setprecision(precision)<<t4<<":"<<t4/iter_time_woutprint*100<<endl;
                    cout<<"save time: "<<setprecision(precision)<<t5<<":"<<t5/iter_time_woutprint*100<<endl;
                    cout<<"print time: "<<setprecision(precision)<<t6<<endl;
                }
            }

            // Write log files
            mainlog(p);
            maxlog(p);
            solverlog(p);
        }

        p->gctime=0.0;
        p->xtime=0.0;
        p->reinitime=0.0;
        p->wavecalctime=0.0;

        #if USE_AMREX
        a->press.FillBoundary();
        #else
        pgc->gcparax(p,a->press,4);
        #endif

        stop(p,a,pgc);
    }

    if(p->mpirank==0)
    {
        cout<<endl<<"******************************"<<endl<<endl;

        cout<<"modelled time: "<<p->simtime<<endl;
        cout << endl;

        mainlogout.close();
        maxlogout.close();
        solvlogout.close();
    }

    pgc->final();
}
