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

#include"initialize.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside_ls.h"
#include "definitions_amrex.h"
#include<utility>

void initialize::iniphi(lexer* p, fdm* a, ghostcell* pgc)
{
    a->phi.setVal(-1.0);

    pgc->start4(p,a->phi,50);

    if(p->F50_flag==1)
    LOOP
    if(p->XN[IP]>=p->F51 && p->XN[IP]<p->F54 && p->YN[JP]>=p->F52 && p->YN[JP]<p->F55 && p->ZN[KP]>=p->F53 && p->ZN[KP]<p->F56)
        a->phi(i,j,k)=1.0;


    if(p->F57_1>0 || p->F57_2>0 || p->F57_3>0 || p->F57_4>0)
    {
        LOOP
        if(p->F57_1*p->XP[IP]+ p->F57_2*p->YP[JP]+ p->F57_3*p->ZP[KP] < p->F57_4)
            a->phi(i,j,k)=1.0;
    }

    if(p->F58_4>0.0)
    {
        double r;
        LOOP
        {
            r = sqrt(pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));

            if(r<=p->F58_4)
                a->phi(i,j,k)=1.0;
        }
    }

    if(p->F59_r>0.0)
    {
        double r;
        LOOP
        {
            r = sqrt(pow(p->XP[IP]-p->F59_xm,2.0)+pow(p->YP[JP]-p->F59_ym,2.0));

            if(r<=p->F59_r && p->pos_z()>p->F59_zs && p->pos_z()<=p->F59_ze)
                a->phi(i,j,k)=1.0;
        }
    }

    if(p->F60>-1.0e20)
    {
        LOOP
            a->phi(i,j,k)=p->F60-p->pos_z();

    }

    if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20 && p->F63>-1.0e-20 )
    {
        double phidiff=p->F62-p->phimean;
        double xdiff=p->xcoormax-p->F63;

        LOOP
        if(p->pos_x() > p->F63)
            a->phi(i,j,k)=(phidiff/xdiff)*(p->pos_x()-p->F63) + p->phimean - p->pos_z() ;
    }

    iniphi_box(p,a,pgc);

    iniphi_wedge(p,a,pgc);

    pgc->start4(p,a->phi,50);

    if(p->F60 > -1.0e20)
    {
        p->phiout=p->F60;
        p->fsfin=p->F60;
        p->fsfout=p->F60;
    }

    if(p->F61 > -1.0e20)
    {
        p->phiin=p->F61;
        p->fsfin=p->F61;
    }

    if(p->F62 > -1.0e20)
    {
        p->phiout=p->F62;
        p->fsfout=p->F62;
    }

    // Initialize ro and visc based on phi
    double H=0.0;
    BASELOOP
    {
        H = heaviside_ls(a->phi(i,j,k), p->psi);

        a->ro(i,j,k)=p->W1*H + p->W3*(1.0-H);
        a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
    }

    pgc->start4(p,a->ro,1);
    pgc->start4(p,a->visc,1);

    // p->level=0;
    // if(p->mpirank==0)std::cout<<"psi: "<<p->psi<<std::endl;
    // for (amrex::MFIter _tile_mfi(p->amr_cell_mf[p->level]/*,MFIter_TILING*/); _tile_mfi.isValid(); ++_tile_mfi) for (struct { lexer* ctx; amrex::MFIter* saved; } _guard{p, p->set_tile_mfi(&_tile_mfi)}; _guard.ctx != nullptr; _guard.ctx->set_tile_mfi(_guard.saved ? _guard.saved : _guard.ctx->default_cell_mfi.get()), _guard.ctx = nullptr)
    // IJKLOOP
    // {
    //     if(k==8&&i==4&&p->mpirank==0)
    //     {
    //         std::cout<<"L0: j: "<<j<<" rho: "<<std::setprecision(6)<<a->ro(i,j,k)<<" phi: "<<a->phi(i,j,k)<<" press: "<<a->press(i,j,k)<<std::endl;
    //     }
    // }
    // p->level=1;
    // int count=0;
    // double temp=0.0, temp2=0.0, temp3=0.0;
    // for (amrex::MFIter _tile_mfi(p->amr_cell_mf[p->level]/*,MFIter_TILING*/); _tile_mfi.isValid(); ++_tile_mfi) for (struct { lexer* ctx; amrex::MFIter* saved; } _guard{p, p->set_tile_mfi(&_tile_mfi)}; _guard.ctx != nullptr; _guard.ctx->set_tile_mfi(_guard.saved ? _guard.saved : _guard.ctx->default_cell_mfi.get()), _guard.ctx = nullptr)
    // IJKLOOP
    // {
    //     if(k+p->amr_tile_lo.z>=16&&k+p->amr_tile_lo.z<=17&&i+p->amr_tile_lo.x==8)
    //     {
    //         temp+=a->ro(i,j,k);
    //         temp2+=a->phi(i,j,k);
    //         temp3+=a->press(i,j,k);
    //         count++;
    //         if(count==4)
    //         {
    //             std::cout<<"L1: j: "<<4+int(j/2)<<" rho average: "<<std::setprecision(6)<<temp/4.0<<" phi: "<<temp2/4.0<<" press: "<<temp3/4.0<<std::endl;
    //             count=0;
    //             temp=0.0;
    //             temp2=0.0;
    //             temp3=0.0;
    //         }
    //     }
    // }
    // p->level=0;
}

void initialize::iniphi_io(fdm*a, lexer* p, ghostcell* pgc)
{
    if(p->F61>-1.0e20)
    GC4LOOP
    {
        GCB4_TILE(n);

        if(p->gcb4[p->level][n].bc==1)
        {
            i=p->gcb4[p->level][n].i;
            j=p->gcb4[p->level][n].j;
            k=p->gcb4[p->level][n].k;

            a->phi(i-1,j,k)=p->F61-p->pos_z();
            a->phi(i-2,j,k)=p->F61-p->pos_z();
            a->phi(i-3,j,k)=p->F61-p->pos_z();
        }
        p->phiin=p->F61;
    }
    GC_TILE_RESET;

    /*if(p->F62>-1.0e20)
    GC4LOOP
    {
        GCB4_TILE(n);

        if(p->gcb4[p->level][n].bc==2)
        {
            i=p->gcb4[p->level][n].i;
            j=p->gcb4[p->level][n].j;
            k=p->gcb4[p->level][n].k;

            a->phi(i+1,j,k)=p->F62-p->pos_z();
            a->phi(i+2,j,k)=p->F62-p->pos_z();
            a->phi(i+3,j,k)=p->F62-p->pos_z();
        }
    }*/
    GC_TILE_RESET;
}

void initialize::iniphi_box(lexer* p, fdm *a, ghostcell* pgc)
{
    int istart, iend, jstart, jend, kstart, kend;
    int qn;

    if(p->F70>0)
    for(qn=0;qn<p->F70;++qn)
    {
        istart = p->posc_i(p->F70_xs[qn]);
        iend = p->posc_i(p->F70_xe[qn]);

        jstart = p->posc_j(p->F70_ys[qn]);
        jend = p->posc_j(p->F70_ye[qn]);

        kstart = p->posc_k(p->F70_zs[qn]);
        kend = p->posc_k(p->F70_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
            a->phi(i,j,k)=1;
    }

    if(p->F71>0)
    for(qn=0;qn<p->F71;++qn)
    {
        istart = p->posc_i(p->F71_xs[qn]);
        iend = p->posc_i(p->F71_xe[qn]);

        jstart = p->posc_j(p->F71_ys[qn]);
        jend = p->posc_j(p->F71_ye[qn]);

        kstart = p->posc_k(p->F71_zs[qn]);
        kend = p->posc_k(p->F71_ze[qn]);


        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
            a->phi(i,j,k)=-1;
    }

    if(p->F72>0)
    for(qn=0;qn<p->F72;++qn)
    {
        istart = p->posc_i(p->F72_xs[qn]);
        iend = p->posc_i(p->F72_xe[qn]);

        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

        LOOP
        if(i>=istart && i<iend && j>=jstart && j<jend)
            a->phi(i,j,k)=p->F72_h[qn]-p->pos_z();
    }
}

void initialize::iniphi_wedge(lexer* p, fdm *a, ghostcell* pgc)
{
    if(p->F112>0)
    {
        for(int qn=0; qn<p->F112; ++qn)
        {
            double slope=(p->F112_ze[qn]-p->F112_zs[qn])/(p->F112_xe[qn]-p->F112_xs[qn]);
            double z = p->F112_zs[qn];
            if(p->F112_ze[qn]<p->F112_zs[qn])
            {
                std::swap(p->F112_ze[qn],p->F112_zs[qn]);
                z = p->F112_ze[qn];
            }

            LOOP
                if(p->pos_x()>=p->F112_xs[qn] && p->pos_x()<p->F112_xe[qn]
                    && p->pos_y()>=p->F112_ys[qn] && p->pos_y()<p->F112_ye[qn]
                    && p->pos_z()>=p->F112_zs[qn] && p->pos_z()<slope*(p->pos_x()-p->F112_xs[qn])+z)
                {
                    a->phi(i,j,k)=1.0;
                }
        }
    }

    if(p->F113>0)
    {
        for(int qn=0; qn<p->F113; ++qn)
        {
            double slope=(p->F113_ze[qn]-p->F113_zs[qn])/(p->F113_ye[qn]-p->F113_ys[qn]);
            double z = p->F113_zs[qn];
            if(p->F113_ze[qn]<p->F113_zs[qn])
            {
                std::swap(p->F113_ze[qn],p->F113_zs[qn]);
                z = p->F113_ze[qn];
            }

            LOOP
                if(p->pos_x()>=p->F113_xs[qn] && p->pos_x()<p->F113_xe[qn]
                    && p->pos_y()>=p->F113_ys[qn] && p->pos_y()<p->F113_ye[qn]
                    && p->pos_z()>=p->F113_zs[qn] && p->pos_z()<slope*(p->pos_y()-p->F113_ys[qn])+z)
                {
                    a->phi(i,j,k)=1.0;
                }
        }
    }

    if(p->F114>0)
    {
        for(int qn=0; qn<p->F114; ++qn)
        {
            double slope=(p->F114_ze[qn]-p->F114_zs[qn])/(p->F114_xe[qn]-p->F114_xs[qn]);
            double z = p->F114_zs[qn];
            if(p->F114_ze[qn]<p->F114_zs[qn])
            {
                std::swap(p->F114_ze[qn],p->F114_zs[qn]);
                z = p->F114_ze[qn];
            }

            LOOP
                if(p->pos_x()>=p->F114_xs[qn] && p->pos_x()<p->F114_xe[qn]
                    && p->pos_y()>=p->F114_ys[qn] && p->pos_y()<p->F114_ye[qn]
                    && p->pos_z()>=z+slope*(p->pos_x()-p->F114_xs[qn]) && p->pos_z()<p->F114_ze[qn])
                {
                    a->phi(i,j,k)=1.0;
                }
        }
    }

    if(p->F115>0)
    {
        for(int qn=0; qn<p->F115; ++qn)
        {
            double slope=(p->F115_ze[qn]-p->F115_zs[qn])/(p->F115_ye[qn]-p->F115_ys[qn]);
            double z = p->F115_zs[qn];
            if(p->F115_ze[qn]<p->F115_zs[qn])
            {
                std::swap(p->F115_ze[qn],p->F115_zs[qn]);
                z = p->F115_ze[qn];
            }

            LOOP
                if(p->pos_x()>=p->F115_xs[qn] && p->pos_x()<p->F115_xe[qn]
                    && p->pos_y()>=p->F115_ys[qn] && p->pos_y()<p->F115_ye[qn]
                    && p->pos_z()>=z+slope*(p->pos_y()-p->F115_ys[qn]) && p->pos_z()<p->F115_ze[qn])
                {
                    a->phi(i,j,k)=1.0;
                }
        }
    }
}

void initialize::iniphi_surfarea(lexer* p, fdm *a, ghostcell* pgc)
{
    double dx,dy,dz,dnorm,dirac;
    double area=0.0;
    double epsi = 1.6*p->DXM;

    LOOP
    {
        epsi = (1.6/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        dx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
        dy = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
        dz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);

        dnorm = sqrt(p->DXN[IP]*p->DXN[IP] + p->DYN[JP]*p->DYN[JP] + p->DZN[KP]*p->DZN[KP]);


        dirac=0.0;

        if(fabs(a->phi(i,j,k))<epsi)
            dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));

        area +=  pow(p->DXM,3.0) * dirac *dnorm;
    }

    area = pgc->globalsum(area);
}
