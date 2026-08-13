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

#include"ioflow_gravity.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"vrans_nhflow.h"
#include"patchBC_interface.h"

ioflow_gravity::ioflow_gravity(lexer *p, ghostcell *pgc, patchBC_interface *ppBC)
{
    pBC = ppBC;

    omega_x = 2.0*PI*p->B191_2;
    omega_y = 2.0*PI*p->B192_2;

    theta_x = p->B191_1*(PI/180.0);
    theta_y = p->B192_1*(PI/180.0);
}

void ioflow_gravity::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    #if USE_AMREX
    const int nlevs = p->nlevs;
    #else
    const int nlevs = 1;
    #endif

    p->gcin.resize_levels(nlevs);
    p->gcout.resize_levels(nlevs);
    int cs = 0;

    LEVEL_LOOP
    GC4LOOP
    {
        GCB4_TILE(n);

        i = p->gcb4[p->level][n].i;
        j = p->gcb4[p->level][n].j;
        k = p->gcb4[p->level][n].k;
        cs = p->gcb4[p->level][n].cs;

        if(p->gcb4[p->level][n].bc==1)
        p->gcin[p->level].push_back({i, j, k, cs});
        else if(p->gcb4[p->level][n].bc==2)
        p->gcout[p->level].push_back({i, j, k, cs});
    }
    GC_TILE_RESET;

    p->gcin_count=p->gcin[0].size();
    p->gcout_count=p->gcout[0].size();
}

void ioflow_gravity::discharge(lexer *p, fdm* a, ghostcell* pgc)
{
    // patchBC
    pBC->patchBC_discharge(p,a,pgc);
}

void ioflow_gravity::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    double omega;
    a->gi = p->W20;
    a->gj = p->W21;
    a->gk = p->W22;


    // -------
    // translation
    if(p->B181==1)
    a->gi += p->B181_1*pow(2.0*PI*p->B181_2,2.0)*sin((2.0*PI*p->B181_2)*p->simtime + p->B181_3*(PI/180.0));

    if(p->B182==1)
    a->gj += p->B182_1*pow(2.0*PI*p->B182_2,2.0)*sin((2.0*PI*p->B182_2)*p->simtime + p->B182_3*(PI/180.0));

    if(p->B183==1)
    a->gk += p->B183_1*pow(2.0*PI*p->B183_2,2.0)*sin((2.0*PI*p->B183_2)*p->simtime + p->B183_3*(PI/180.0));


    // -------
    // rotation
    if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        a->gj += sin(theta_x*sin(omega_x*p->simtime))*p->W22;

        a->gk +=  cos(theta_x*sin(omega_x*p->simtime))*p->W22 - p->W22;
    }

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        a->gi += sin(theta_y*sin(omega_y*p->simtime))*p->W22;

        a->gk +=  cos(theta_y*sin(omega_y*p->simtime))*p->W22  - p->W22;
    }

    // -------
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_gravity::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_gravity::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_gravity::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    pBC->patchBC_waterlevel(p,a,pgc,f);
}

void  ioflow_gravity::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    NLOOP4
    a->rhsvec.V[n]=0.0;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        ULOOP
        {
            dist_x = p->pos_x() - p->B192_3;
            dist_z = p->pos_z() - p->B192_4;
            a->rhsvec.V[n] += dist_z*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
                         + dist_x*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
                         //- a->w(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);
            ++n;
        }

        //a->gi = sin(theta_y*sin(omega_y*p->simtime))*p->W22;
        /*a->gi = p->W22*sin( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
                            + pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
                                + 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

    }
}

void  ioflow_gravity::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    NLOOP4
    a->rhsvec.V[n]=0.0;

    //if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    //a->gj = sin(theta_x*sin(omega_x*p->simtime))*p->W22;
}

void  ioflow_gravity::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    NLOOP4
    a->rhsvec.V[n]=0.0;

    if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        WLOOP
        {
            dist_y = p->pos_y() - p->B191_3;
            a->rhsvec.V[n] -= dist_y*theta_x*pow(omega_x,2.0)*sin(omega_x*p->simtime);

            ++n;
        }

        //a->gk =  cos(theta_x*sin(omega_x*p->simtime))*p->W22;
    }

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        WLOOP
        {
            dist_x = p->pos_x() - p->B192_3;
            dist_z = p->pos_z() - p->B192_4;
            a->rhsvec.V[n] +=  -dist_x*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
                         +  dist_z*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
                         //- a->u(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);
            ++n;
        }

        //a->gk =  cos(theta_y*sin(omega_y*p->simtime))*p->W22;
        /*
        a->gk = p->W22*cos( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
                        + pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
                            + 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

    }
}

void  ioflow_gravity::isource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans_nhflow *pvrans, slice &WL)
{
    NLOOP4
    d->rhsvec.V[n]=0.0;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        LOOP
        {
            dist_x = p->pos_x() - p->B192_3;
            dist_z = p->pos_z() - p->B192_4;
            d->rhsvec.V[n] += dist_z*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
                         + dist_x*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
            ++n;
        }

        d->gi = sin(theta_y*sin(omega_y*p->simtime))*p->W22;
    }
    
    //VRANS
    pvrans->u_source(p,d,WL);
}

void  ioflow_gravity::jsource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans_nhflow *pvrans, slice &WL)
{
    NLOOP4
    d->rhsvec.V[n]=0.0;

    if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    d->gj = sin(theta_x*sin(omega_x*p->simtime))*p->W22;
    
    //VRANS
    pvrans->v_source(p,d,WL);
}

void  ioflow_gravity::ksource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans_nhflow *pvrans, slice &WL)
{
    NLOOP4
    d->rhsvec.V[n]=0.0;

    if(p->B191==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        LOOP
        {
            dist_y = p->pos_y() - p->B191_3;
            d->rhsvec.V[n] -= dist_y*theta_x*pow(omega_x,2.0)*sin(omega_x*p->simtime);

            ++n;
        }

        //d->gk =  cos(theta_x*sin(omega_x*p->simtime))*p->W22;
    }

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    {
        n=0;
        WLOOP
        {
            dist_x = p->pos_x() - p->B192_3;
            dist_z = p->pos_z() - p->B192_4;
            d->rhsvec.V[n] +=  -dist_x*theta_y*pow(omega_y,2.0)*sin(omega_y*p->simtime)
                         +  dist_z*pow(theta_y*omega_y*cos(omega_y*p->simtime),2.0);
                         //- a->u(i,j,k)*theta_y*omega_y*cos(omega_y*p->simtime);
            ++n;
        }

        //d->gk =  cos(theta_y*sin(omega_y*p->simtime))*p->W22;
        /*
        a->gk = p->W22*cos( theta_y*sin(omega_y*p->simtime) - 0.31*theta_y*theta_y*(1.0+cos(2.0*omega_y*p->simtime))
                        + pow(theta_y,3.0)*(0.16*cos(omega_y*p->simtime) - 0.16*cos(3.0*omega_y*p->simtime)
                            + 0.13*sin(omega_y*p->simtime) - 0.004*sin(3.0*omega_y*p->simtime)));*/

    }
    
    //VRANS
    pvrans->w_source(p,d,WL);
}


void ioflow_gravity::pressure_io(lexer *p, fdm *a, ghostcell *pgc)
{
    pBC->patchBC_pressure(p,a,pgc,a->press);
}

double ioflow_gravity::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    return 0.0;
}

double ioflow_gravity::wave_xvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_gravity::wave_yvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_gravity::wave_zvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

int ioflow_gravity::iozonecheck(lexer *p, fdm*a)
{
    return 1;
}

void ioflow_gravity::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_gravity::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_gravity::isource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    SLICELOOP1
    b->F(i,j)=0.0;
}

void ioflow_gravity::jsource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    SLICELOOP2
    b->G(i,j)=0.0;
}

void ioflow_gravity::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    ALOOP
    a->porosity(i,j,k)=1.0;

    pgc->start4(p,a->porosity,1);
}
