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

#include"ioflow_void.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"vrans.h"
#include"rheology_v.h"
#include"rheology_f.h"
#include"turbulence.h"
#include"patchBC_interface.h"

ioflow_v::ioflow_v(lexer *p, ghostcell *pgc, patchBC_interface *ppBC) : flowfile_in(p,pgc)
{
    pBC = ppBC;
}

ioflow_v::~ioflow_v()
{
    delete prheo;
}

void ioflow_v::discharge(lexer *p, fdm* a, ghostcell* pgc)
{
    // patchBC
    pBC->patchBC_discharge(p,a,pgc);
}

void ioflow_v::inflow(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    if(p->I230>0)
    ff_inflow(p,a,pgc,u,v,w);

    prheo->filltau(p,a,pgc);

    velocity_inlet(p,a,pgc,u,v,w);

    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_v::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    velocity_inlet(p,a,pgc,u,v,w);

    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void ioflow_v::velocity_inlet(lexer *p, fdm* a, ghostcell* pgc, field &u, field &v, field &w)
{
    if(p->W11==1 && p->bc_type[0] == 1)
    {
        u.FillDomainBoundaryValue(p->W11_u, 0, false);
        v.FillDomainBoundaryValue(p->W11_v, 0, false);
        w.FillDomainBoundaryValue(p->W11_w, 0, false);
    }
    if(p->W14==1 && p->bc_type[1] == 1)
    {
        u.FillDomainBoundaryValue(p->W14_u, 0, true);
        v.FillDomainBoundaryValue(p->W14_v, 0, true);
        w.FillDomainBoundaryValue(p->W14_w, 0, true);
    }
    if(p->W13==3 && p->bc_type[2] == 1)
    {
        u.FillDomainBoundaryValue(p->W13_u, 1, false);
        v.FillDomainBoundaryValue(p->W13_v, 1, false);
        w.FillDomainBoundaryValue(p->W13_w, 1, false);
    }
    if(p->W12==1 && p->bc_type[3] == 1)
    {
        u.FillDomainBoundaryValue(p->W12_u, 1, true);
        v.FillDomainBoundaryValue(p->W12_v, 1, true);
        w.FillDomainBoundaryValue(p->W12_w, 1, true);
    }
    if(p->W15==1 && p->bc_type[4] == 1)
    {
        u.FillDomainBoundaryValue(p->W15_u, 2, false);
        v.FillDomainBoundaryValue(p->W15_v, 2, false);
        w.FillDomainBoundaryValue(p->W15_w, 2, false);
    }
    if(p->W16==1 && p->bc_type[5] == 1)
    {
        u.FillDomainBoundaryValue(p->W16_u, 2, true);
        v.FillDomainBoundaryValue(p->W16_v, 2, true);
        w.FillDomainBoundaryValue(p->W16_w, 2, true);
    }
}

void ioflow_v::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->I230>0)
    ff_waterlevel(p,a,pgc,a->phi);

    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_v::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    pBC->patchBC_waterlevel(p,a,pgc,f);
}

void ioflow_v::iogcb_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void  ioflow_v::isource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    a->rhsvec.V[n]=0.0;

    count=0;
    if(p->B240>0 && p->B241==1)
    ULOOP
    {
        // porous media
        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*a->visc(i,j,k)*a->u(i,j,k) + 0.5*p->B240_C[n]*a->u(i,j,k)*fabs(a->u(i,j,k));
        }

        a->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    pvrans->u_source(p,a);

    //Rheology
    prheo->u_source(p,a);
}

void  ioflow_v::jsource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    a->rhsvec.V[n]=0.0;

    count=0;
    if(p->B240>0 && p->B242==1)
    VLOOP
    {
        // porous media
        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*a->visc(i,j,k)*a->v(i,j,k) + 0.5*p->B240_C[n]*a->v(i,j,k)*fabs(a->v(i,j,k));
        }

        a->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    pvrans->v_source(p,a);

    //Rheology
    prheo->v_source(p,a);
}

void  ioflow_v::ksource(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    a->rhsvec.V[n]=0.0;

    count=0;
    if(p->B240>0 && p->B243==1)
    WLOOP
    {
        // porous media
        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*a->visc(i,j,k)*a->w(i,j,k) + 0.5*p->B240_C[n]*a->w(i,j,k)*fabs(a->w(i,j,k));
        }

        a->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    pvrans->w_source(p,a);

    //Rheology
    prheo->w_source(p,a);
}

void ioflow_v::isource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    d->rhsvec.V[n]=0.0;

    // Darcy Porosity
    count=0;
    if(p->B240>0 && p->B241==1)
    LOOP
    {

        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*d->VISC[IJK]*d->U[IJK] + 0.5*p->B240_C[n]*d->U[IJK]*fabs(d->U[IJK]);
        }

        d->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    //pvrans->u_source(p,a);
}

void ioflow_v::jsource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    d->rhsvec.V[n]=0.0;

    count=0;
    if(p->B240>0 && p->B242==1)
    VLOOP
    {
        // porous media
        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*d->VISC[IJK]*d->V[IJK] + 0.5*p->B240_C[n]*d->V[IJK]*fabs(d->V[IJK]);
        }

        d->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    //pvrans->v_source(p,a);
}

void ioflow_v::ksource_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, vrans *pvrans)
{
    double porousterm;

    NLOOP4
    d->rhsvec.V[n]=0.0;

    count=0;
    if(p->B240>0 && p->B243==1)
    LOOP
    {
        // porous media
        porousterm=0.0;
        for(n=0;n<p->B240;++n)
        {
            if(p->pos_x() >= p->B240_xs[n] && p->pos_x() < p->B240_xe[n])
            if(p->pos_y() >= p->B240_ys[n] && p->pos_y() < p->B240_ye[n])
            if(p->pos_z() >= p->B240_zs[n] && p->pos_z() < p->B240_ze[n])
            porousterm=p->B240_D[n]*d->VISC[IJK]*d->W[IJK] + 0.5*p->B240_C[n]*d->W[IJK]*fabs(d->W[IJK]);
        }

        d->rhsvec.V[count] -= porousterm;
        ++count;
    }

    //VRANS
    //pvrans->w_source(p,a);
}

void ioflow_v::pressure_io(lexer *p, fdm *a, ghostcell* pgc)
{
    double pval=0.0;

    GC4LOOP
    if(p->gcb4[n][4]==2)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        pval=0.0;

        if(p->B77==1)
        {
            pval=(p->phiout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
        }

        if(p->B77==10)
        {
            double eps,H;

            eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

            if(a->phi(i,j,k)>eps)
            H=1.0;
            if(a->phi(i,j,k)<-eps)
            H=0.0;
            if(fabs(a->phi(i,j,k))<=eps)
            H=0.5*(1.0 + a->phi(i,j,k)/eps + (1.0/PI)*sin((PI*a->phi(i,j,k))/eps));

            pval=(1.0-H)*a->press(i,j,k);

            a->press(i,j,k)=pval;
            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
        }
    }

    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void ioflow_v::u_relax(lexer *p, fdm *a, ghostcell *pgc, field &uvel)
{
    double epsi,H,fbval;
    double dist;
    double cosb,sinb;

    for(int qn=0; qn<p->W41; ++qn)
    {
        cosb = cos(p->W41_beta[qn]*PI/180.0);

        ULOOP
        {
            dist = sqrt(pow(p->W41_xc[qn]-p->pos1_x(),2.0) + pow(p->W41_yc[qn]-p->pos_y(),2.0));

            if(dist>epsi)
            H=1.0;
            if(dist<-epsi)
            H=0.0;
            if(fabs(dist)<=epsi)
            H=0.5*(1.0 + dist/epsi + (1.0/PI)*sin((PI*dist)/epsi));

            if(0.5*(a->phi(i,j,k)+a->phi(i+1,j,k))>0.0)
            a->u(i,j,k) = H*a->u(i,j,k) + (1.0-H)*p->W41_vel[qn]*cosb;
        }
    }

}

void ioflow_v::v_relax(lexer *p, fdm *a, ghostcell *pgc, field &vvel)
{
    double epsi,H,fbval;
    double dist;
    double cosb,sinb;

    for(int qn=0; qn<p->W41; ++qn)
    {
        sinb = sin(p->W41_beta[qn]*PI/180.0);

        VLOOP
        {
            dist = sqrt(pow(p->W41_xc[qn]-p->pos_x(),2.0) + pow(p->W41_yc[qn]-p->pos2_y(),2.0));

            if(dist>epsi)
            H=1.0;
            if(dist<-epsi)
            H=0.0;
            if(fabs(dist)<=epsi)
            H=0.5*(1.0 + dist/epsi + (1.0/PI)*sin((PI*dist)/epsi));

            if(0.5*(a->phi(i,j,k)+a->phi(i,j+1,k))>0.0)
            a->v(i,j,k) = H*a->v(i,j,k) + (1.0-H)*p->W41_vel[qn]*sinb;
        }
    }
}

double ioflow_v::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    return 0.0;
}

double ioflow_v::wave_xvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_v::wave_yvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_v::wave_zvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

int ioflow_v::iozonecheck(lexer *p, fdm*a)
{
    return 1;
}

void ioflow_v::discharge2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    pBC->patchBC_discharge2D(p,b,pgc,b->P,b->Q,b->eta,b->bed);
}

void ioflow_v::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_v::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_v::isource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    SLICELOOP1
    b->F(i,j)=0.0;
}

void ioflow_v::jsource2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    SLICELOOP2
    b->G(i,j)=0.0;
}

void ioflow_v::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    if(p->W90>0)
    prheo = new rheology_f(p);
    else
    prheo = new rheology_v();
}

void ioflow_v::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc,vrans *pvrans)
{
    pvrans->sed_update(p,a,pgc);
}

void ioflow_v::waterlevel2D(lexer *p, fdm2D *b, ghostcell* pgc, slice &eta)
{
    pBC->patchBC_waterlevel2D(p,b,pgc,eta);
}
