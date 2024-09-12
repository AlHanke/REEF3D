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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include<math.h>

void partres::move_RK2_step1(lexer *p, fdm &a, ghostcell &pgc, particles_obj2 &PP, sediment_fdm &s, turbulence &pturb, int &xchanged, int &removed)
{
    particles_obj2 Send[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};
    particles_obj2 Recv[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};
    
    particlePerCell(p,pgc,PP);
    particleStressTensor(p,a,pgc,PP);
    timestep(p,pgc,PP);
    
    Umax=-1.0e10;
    
    // RK step 1 
    for(auto &particle : PP.particles)
    if(particle.flag>=0)
    {
        if(p->Q11==1)
        advec_plain(p, a, particle, PP.d50, s, pturb, 
                        particle.x, particle.y, particle.z, particle.u, particle.v, particle.w,
                        F, G, H, 1.0);
        
        // if(p->Q11==2)
        // advec_pic(p, a, PP, n, s, pturb, 
        //                 particle.X, particle.Y, particle.Z, particle.U, particle.V, particle.W,
        //                 F, G, H, 1.0);
                                         
        // Velocity update
        particle.urk1 = particle.u + p->dtsed*F;
        particle.vrk1 = particle.v + p->dtsed*G;
        particle.wrk1 = 0.0; //particle.w + p->dtsed*H;
        
        // Position update
        particle.xrk1 = particle.x + p->dtsed*particle.urk1;
        particle.yrk1 = particle.y + p->dtsed*particle.vrk1;
        particle.zrk1 = particle.z + p->dtsed*particle.wrk1;
        
        //if(particle.Flag[n]==1)
        //particle.zrk1 = p->ccslipol4(s.bedzh,particle.xrk1,particle.yrk1);

        // Particel sum update
        cellSum[IJK]-=particle.parcelFactor;
        bedChange[IJ]-=particle.parcelFactor;
        i=p->posc_i(particle.xrk1);
        j=p->posc_j(particle.yrk1);
        k=p->posc_k(particle.zrk1);
        cellSum[IJK]+=particle.parcelFactor;
        bedChange[IJ]+=particle.parcelFactor;
        particleStressTensorUpdateIJK(p,a,PP);

        addParticleForTransfer(p,PP,n,Send,xchanged);
    }
    
    // vertical coordinate
    Umax = pgc.globalmax(Umax);
    
    for(auto &particle : PP.particles)
    if(particle.flag>=0)
    {
        Uabs=sqrt(particle.urk1*particle.urk1+particle.vrk1*particle.vrk1);
        
        fac = Uabs/(Umax>1.0e-10?Umax:1.0e10);
        
        //cout<<"Uabs: "<<Uabs<<" Umax: "<<Umax<<" fac: "<<fac<<endl;
        
        if(Uabs>=0.1*Umax && Uabs>0.01)
        {
        k=p->posc_k(particle.zrk1);
        particle.zrk1 =   s.bedzh(i,j) + 0.5*p->DZN[KP]*double(rand() % irand)/drand;
            
        }
    }

    // recv
    {
        pgc.para_particles_obj2(p,Send,Recv);

        for(int n=0;n<4;n++)
        {
            for(auto &particle : Recv[n].particles)
            {
                i = p->posc_i(particle.xrk1);
                j = p->posc_j(particle.yrk1);
                k = p->posc_k(particle.zrk1);
                transfer(p,particle);
            }
            PP.add_obj(&Recv[n]);
        }
    }


    {
        bool inBounds=false;
        int i,j,k;
        boundarycheck bounderies;

        for(auto &particle : PP.particles)
            if(particle.flag>=0)
            {
                i = p->posc_i(particle.xrk1);
                j = p->posc_j(particle.yrk1);
                k = p->posc_k(particle.zrk1);

                inBounds=bounderies.globalminboundcheck(p,i,j,k);
                if (inBounds)
                    inBounds=bounderies.globalmaxboundcheck(p,i,j,k);

                // remove out of bounds particles
                if(!inBounds)
                {
                    //remove(p,PP,n);
                    PP.erase(n);
                    removed++;
                }
            }
    }
    
    particleStressTensor(p,a,pgc,PP);
}
    
void partres::move_RK2_step2(lexer *p, fdm &a, ghostcell &pgc, particles_obj2 &PP, sediment_fdm &s, turbulence &pturb, int &xchanged, int &removed)
{
    particles_obj2 Send[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};
    particles_obj2 Recv[6]={particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),
    particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density),particles_obj2(10,PP.d50,PP.density)};
    
    Umax=-1.0e10;
    
    // RK step 2
    for(auto &particle : PP.particles)
    if(particle.flag>=0)
    {
        if(p->Q11==1)
        advec_plain(p, a, particle, PP.d50, s, pturb, 
            particle.xrk1, particle.yrk1, particle.zrk1, particle.urk1, particle.vrk1, particle.wrk1,
            F, G, H, 0.5);
                        
        // if(p->Q11==2)
        // advec_pic(p, a, PP, n, s, pturb, 
        //                 particle.XRK1, particle.YRK1, particle.ZRK1, particle.URK1, particle.VRK1, particle.WRK1,
        //                 F, G, H, 0.5);
                        
        // Velocity update
        particle.u = 0.5*particle.u + 0.5*particle.urk1 + 0.5*p->dtsed*F;
        particle.v = 0.5*particle.v + 0.5*particle.vrk1 + 0.5*p->dtsed*G;
        particle.w = 0.0;//0.5*particle.w + 0.5*particle.urk1 + 0.5*p->dtsed*H;
        
        // Position update
        particle.x = 0.5*particle.x + 0.5*particle.xrk1 + 0.5*p->dtsed*particle.u;
        particle.y = 0.5*particle.y + 0.5*particle.yrk1 + 0.5*p->dtsed*particle.v;
        particle.z = 0.5*particle.z + 0.5*particle.zrk1 + 0.5*p->dtsed*particle.w;
        
        //if(particle.Flag[n]==1)
        //particle.z = p->ccslipol4(s.bedzh,particle.x,particle.y);

        // Particel sum update
        cellSum[IJK]-=particle.parcelFactor;
        bedChange[IJ]-=particle.parcelFactor;
        i=p->posc_i(particle.x);
        j=p->posc_j(particle.y);
        k=p->posc_k(particle.z);
        cellSum[IJK]+=particle.parcelFactor;
        bedChange[IJ]+=particle.parcelFactor;
        particleStressTensorUpdateIJK(p,a,PP);

        addParticleForTransfer(p,PP,n,Send,xchanged);
    }
    
    // vertical coordinate
    Umax = pgc.globalmax(Umax);
    
    cout<<"Umax Particle: "<<Umax<<endl;
    
    for(auto &particle : PP.particles)
    if(particle.flag>=0)
    {
        Uabs=sqrt(particle.u*particle.u+particle.v*particle.v);
        
        fac = Uabs/(Umax>1.0e-10?Umax:1.0e10);
    
        if(Uabs>=0.1*Umax && Uabs>0.01)
        {
        k=p->posc_k(particle.z);
        particle.z =   s.bedzh(i,j) + 0.5*p->DZN[KP]*double(rand() % irand)/drand;
            
        }
    }
    
    
    {
        pgc.para_particles_obj2(p,Send,Recv);

        for(int n=0;n<4;n++)
        {
            for(auto &particle : Recv[n].particles)
            {
                i = p->posc_i(particle.x);
                j = p->posc_j(particle.y);
                k = p->posc_k(particle.z);
                transfer(p,particle);
            }
            PP.add_obj(&Recv[n]);
        }
    }


    {
        bool inBounds=false;
        int removed=0;
        int i,j,k;
        boundarycheck bounderies;

        for(auto &particle : PP.particles)
            if(particle.flag>=0)
            {
                i = p->posc_i(particle.x);
                j = p->posc_j(particle.y);
                k = p->posc_k(particle.z);

                inBounds=bounderies.globalminboundcheck(p,i,j,k);
                if (inBounds)
                    inBounds=bounderies.globalmaxboundcheck(p,i,j,k);

                // remove out of bounds particles
                if(!inBounds)
                {
                    //remove(p,PP,n);
                    PP.erase(n);
                    removed++;
                }
            }
    }

    if(p->mpirank==0)
    {
        p->sedtime += p->dtsed;
        cout<<"Sediment time: "<<p->sedtime<<" time step: "<<p->dtsed<<endl;
    }
}

void partres::move_RK2(lexer *p, fdm &a, ghostcell &pgc, particles_obj2 &PP, sediment_fdm &s, turbulence &pturb)
{
    // particlePerCell(p,pgc,PP);
    // particleStressTensor(p,a,pgc,PP);
    // timestep(p,pgc,PP);
    
    // // RK step 1 
    // for(size_t n=0;n<particle.loopindex;n++)
    // if(particle.Flag[n]>=0)
    // {
    //     if(p->Q11==1)
    //     advec_plain(p, a, PP, n, s, pturb, 
    //                     particle.X, particle.Y, particle.Z, particle.U, particle.V, particle.W,
    //                     F, G, H, 1.0);
        
    //     if(p->Q11==2)
    //     advec_pic(p, a, PP, n, s, pturb, 
    //                     particle.X, particle.Y, particle.Z, particle.U, particle.V, particle.W,
    //                     F, G, H, 1.0);
                                         
    //     // Velocity update
    //     particle.urk1 = particle.u + p->dtsed*F;
    //     particle.vrk1 = particle.v + p->dtsed*G;
    //     particle.wrk1 = 0.0; //particle.w + p->dtsed*H;
        
    //     // Position update
    //     particle.xrk1 = particle.x + p->dtsed*particle.urk1;
    //     particle.yrk1 = particle.y + p->dtsed*particle.vrk1;
    //     particle.zrk1 = particle.z + p->dtsed*particle.wrk1;

    //     // Particel sum update
    //     cellSum[IJK]-=particle.parcelFactor;
    //     bedChange[IJ]-=particle.parcelFactor;
    //     i=p->posc_i(particle.xrk1);
    //     j=p->posc_j(particle.yrk1);
    //     k=p->posc_k(particle.zrk1);
    //     cellSum[IJK]+=particle.parcelFactor;
    //     bedChange[IJ]+=particle.parcelFactor;
    //     particleStressTensorUpdateIJK(p,a,PP);
    // }
    
    // particleStressTensor(p,a,pgc,PP);
    
    
    // // RK step 2
    // for(size_t n=0;n<particle.loopindex;n++)
    // if(particle.Flag[n]>=0)
    // {
    //     if(p->Q11==1)
    //     advec_plain(p, a, PP, n, s, pturb, 
    //                     particle.XRK1, particle.YRK1, particle.ZRK1, particle.URK1, particle.VRK1, particle.WRK1,
    //                     F, G, H, 0.5);
                        
    //     if(p->Q11==2)
    //     advec_pic(p, a, PP, n, s, pturb, 
    //                     particle.XRK1, particle.YRK1, particle.ZRK1, particle.URK1, particle.VRK1, particle.WRK1,
    //                     F, G, H, 0.5);
                        
    //     // Velocity update
    //     particle.u = 0.5*particle.u + 0.5*particle.urk1 + 0.5*p->dtsed*F;
    //     particle.v = 0.5*particle.v + 0.5*particle.vrk1 + 0.5*p->dtsed*G;
    //     particle.w = 0.0;//0.5*particle.w + 0.5*particle.urk1 + 0.5*p->dtsed*H;
        
    //     // Position update
    //     particle.x = 0.5*particle.x + 0.5*particle.xrk1 + 0.5*p->dtsed*particle.u;
    //     particle.y = 0.5*particle.y + 0.5*particle.yrk1 + 0.5*p->dtsed*particle.v;
    //     particle.z = 0.5*particle.z + 0.5*particle.zrk1 + 0.5*p->dtsed*particle.w;

    //     // Particel sum update
    //     cellSum[IJK]-=particle.parcelFactor;
    //     bedChange[IJ]-=particle.parcelFactor;
    //     i=p->posc_i(particle.x);
    //     j=p->posc_j(particle.y);
    //     k=p->posc_k(particle.z);
    //     cellSum[IJK]+=particle.parcelFactor;
    //     bedChange[IJ]+=particle.parcelFactor;
    //     particleStressTensorUpdateIJK(p,a,PP);
    // }

    // if(p->mpirank==0)
    // {
    //     p->sedtime += p->dtsed;
    //     cout<<"Sediment time: "<<p->sedtime<<" time step: "<<p->dtsed<<endl;
    // }
    
}