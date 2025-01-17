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

#include "particle_func.h"

#include "lexer.h"
#include "fdm.h"

#define PARTICLELOOP for(size_t n=0;n<PP->loopindex;n++) if(PP->Flag[n]>INT32_MIN)

/// @brief Applies advection to positions of particles in @param PP
/// @param p partition object
/// @param a fdm object contains flow field
/// @param PP tracers_obj contains tracer information
/// @param minflag PP.Flag[n] needs to be bigger than minflag for PP[n] to be affected by avection
// void particle_func::advect(lexer* p, fdm* a, tracers_obj* PP, int minflag, double source_u, double source_v, double source_w)
// {
//     double coord1, coord2, coord3, u1, u2, v1, v2, w1, w2;
//     PARTICLELOOP
//         if(PP->Flag[n]>minflag)
//         {
//             i=p->posc_i(PP->X[n]);
//             j=p->posc_j(PP->Y[n]);
//             k=p->posc_k(PP->Z[n]);
//             cellSum[IJK]--;

//             u1=p->dt*(p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n])+source_u);
//             coord1=PP->X[n]+u1;
            
//             v1=p->dt*(p->ccipol2(a->v,PP->X[n],PP->Y[n],PP->Z[n])+source_v);
//             coord2=PP->Y[n]+v1;
            
//             w1=p->dt*(p->ccipol3(a->w,PP->X[n],PP->Y[n],PP->Z[n])+source_w);
//             coord3=PP->Z[n]+w1;
            
            
//             u2=0.25*u1 + 0.25*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);
//             coord1=PP->X[n]+u2;
            
//             v2=0.25*v1 + 0.25*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
//             coord2=PP->Y[n]+v2;
            
//             w2=0.25*w1 + 0.25*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);
//             coord3=PP->Z[n]+w2;
            
            
//             PP->X[n] += (2.0/3.0)*u2 + (2.0/3.0)*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);

//             PP->Y[n] += (2.0/3.0)*v2 + (2.0/3.0)*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
            
//             PP->Z[n] += (2.0/3.0)*w2 + (2.0/3.0)*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);

//             i=p->posc_i(PP->X[n]);
//             j=p->posc_j(PP->Y[n]);
//             k=p->posc_k(PP->Z[n]);
//             cellSum[IJK]++;
//         }
// }

/// @brief Applies advection to positions of particles in @param PP
/// @param p partition object
/// @param a fdm object contains flow field
/// @param PP tracers_obj contains tracer information
/// @param minflag PP.Flag[n] needs to be bigger than minflag for PP[n] to be affected by avection
void particle_func::advect(lexer* p, fdm* a, particles_obj* PP, int minflag, double source_u, double source_v, double source_w)
{
    double coord1, coord2, coord3, u1, u2, v1, v2, w1, w2;
    PARTICLELOOP
        if(PP->Flag[n]>minflag)
        {
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            cellSum[IJK]-=PP->PackingFactor[n];

            source_w -=settling_vel(p,a,PP,n);
            u1=p->dt*(p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n])+source_u);
            coord1=PP->X[n]+u1;
            
            v1=p->dt*(p->ccipol2(a->v,PP->X[n],PP->Y[n],PP->Z[n])+source_v);
            coord2=PP->Y[n]+v1;
            
            w1=p->dt*(p->ccipol3(a->w,PP->X[n],PP->Y[n],PP->Z[n])+source_w);
            coord3=PP->Z[n]+w1;
            
            
            u2=0.25*u1 + 0.25*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);
            coord1=PP->X[n]+u2;
            
            v2=0.25*v1 + 0.25*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
            coord2=PP->Y[n]+v2;
            
            w2=0.25*w1 + 0.25*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);
            coord3=PP->Z[n]+w2;
            
            
            PP->X[n] += (2.0/3.0)*u2 + (2.0/3.0)*p->dt*(p->ccipol1(a->u,coord1,coord2,coord3)+source_u);

            PP->Y[n] += (2.0/3.0)*v2 + (2.0/3.0)*p->dt*(p->ccipol2(a->v,coord1,coord2,coord3)+source_v);
            
            PP->Z[n] += (2.0/3.0)*w2 + (2.0/3.0)*p->dt*(p->ccipol3(a->w,coord1,coord2,coord3)+source_w);

            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            cellSum[IJK]+=PP->PackingFactor[n];
        }
}

/// @brief Particle transport equation following (Tavouktsoglou 2021)
/// @param p 
/// @param a 
/// @param PP 
/// @param flag Particles need to have a larger flag than this
/// Uses Runge Kutta 3 to calculate the change in particle velocity 
/// The particle velocity is then used to transport the particle
// void particle_func::transport(lexer* p, fdm* a, particles_obj* PP, int flag)
// {
//     double RKu,RKv,RKw;
//     double u,v,w;
//     double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
//     /// @brief Difference between flowfield and particle velocity
//     double du, dv, dw;
//     double Dp, thetas;
//     double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
//     double stressDivX=0, stressDivY=0, stressDivZ=0;
//     double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;
//     bool print=true;

//     PARTICLELOOP
//         if(PP->Flag[n]>flag)
//         {
//             // Prep
//             i=p->posc_i(PP->X[n]);
//             j=p->posc_j(PP->Y[n]);
//             k=p->posc_k(PP->Z[n]);

//             thetas=theta_s(p,a,PP,i,j,k);

//             u=p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n]);
//             v=p->ccipol1(a->v,PP->X[n],PP->Y[n],PP->Z[n]);
//             w=p->ccipol1(a->w,PP->X[n],PP->Y[n],PP->Z[n]);

//             stressDivX = (stressTensor[Ip1JK] - stressTensor[IJK])/(p->DXN[IP]);
//             stressDivY = (0.5*(stressTensor[IJp1K]+stressTensor[Ip1Jp1K]) - 0.5*(stressTensor[IJm1K]+stressTensor[Ip1Jm1K]))/(p->DYN[JM1]+p->DYN[JP]);
//             stressDivZ = (0.5*(stressTensor[IJKp1]+stressTensor[Ip1JKp1]) - 0.5*(stressTensor[IJKm1]+stressTensor[Ip1JKm1]))/(p->DYN[KM1]+p->DYN[KP]);

//             // stressDivX = (stressTensor[Ip1JK] - stressTensor[Im1JK])/(p->DXN[IP]+p->DXN[IM1]);
//             stressDivY = (stressTensor[IJp1K] - stressTensor[IJm1K])/(p->DYN[JP]+p->DYN[JM1]);
//             stressDivZ = (stressTensor[IJKp1] - stressTensor[IJKm1])/(p->DZN[KP]+p->DZN[KM1]);

//             pressureDivX = (a->press(i+1,j,k) - a->press(i,j,k))/(p->DXN[IP]);
//             pressureDivY = (0.5*(a->press(i,j+1,k)+a->press(i+1,j+1,k)) - 0.5*(a->press(i,j-1,k)+a->press(i+1,j-1,k)))/(p->DYN[JM1]+p->DYN[JP]);
//             pressureDivZ = (0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1)) - 0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DYN[KM1]+p->DYN[KP]);

//             // RK3 step 1
//             du=u-PP->U[n];
//             dv=v-PP->V[n];
//             dw=w-PP->W[n];

//             Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

//             du1=Dp*du+netBuoyX-pressureDivX/p->S22-stressDivX/((1-thetas)*p->S22);
//             dv1=Dp*dv+netBuoyY-pressureDivY/p->S22-stressDivY/((1-thetas)*p->S22);
//             dw1=Dp*dw+netBuoyZ-pressureDivZ/p->S22-stressDivZ/((1-thetas)*p->S22);

//             RKu=PP->U[n]+du1*p->dt;
//             RKv=PP->V[n]+dv1*p->dt;
//             RKw=PP->W[n]+dw1*p->dt;
            
//             // RK step 2
//             du=u-RKu;
//             dv=v-RKv;
//             dw=w-RKw;

//             Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

//             du2=Dp*du+netBuoyX-(pressureDivX/p->S22+stressDivX/((1-thetas)*p->S22));
//             dv2=Dp*dv+netBuoyY-(pressureDivY/p->S22+stressDivY/((1-thetas)*p->S22));
//             dw2=Dp*dw+netBuoyZ-(pressureDivZ/p->S22+stressDivZ/((1-thetas)*p->S22));

//             du2=0.25*du2+0.25*du1;
//             dv2=0.25*dv2+0.25*dv1;
//             dw2=0.25*dw2+0.25*dw1;

//             RKu=PP->U[n]+du2*p->dt;
//             RKv=PP->V[n]+dv2*p->dt;
//             RKw=PP->W[n]+dw2*p->dt;
            
//             // RK step 3
//             du=u-RKu;
//             dv=v-RKv;
//             dw=w-RKw;

//             Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

//             du3=Dp*du+netBuoyX-(pressureDivX/p->S22+stressDivX/((1-thetas)*p->S22));
//             dv3=Dp*dv+netBuoyY-(pressureDivY/p->S22+stressDivY/((1-thetas)*p->S22));
//             dw3=Dp*dw+netBuoyZ-(pressureDivZ/p->S22+stressDivZ/((1-thetas)*p->S22));


//             if(du2!=du2||du3!=du3)
//             {
//                 cerr<<"Particle velocity component u resulted in NaN.\n"
//                 <<du2<<","<<du3<<"|"<<Dp<<","<<netBuoyX<<","<<pressureDivX<<","<<stressDivX
//                 <<endl;
//                 exit(1);
//             }
//             else
//                 PP->U[n] += ((2.0/3.0)*du2 + (2.0/3.0)*du3)*p->dt;
//             if(dv2!=dv2||dv3!=dv3)
//             {
//                 cerr<<"Particle velocity component v resulted in NaN.\n"
//                 <<dv2<<","<<dv3<<"|"<<Dp<<","<<netBuoyY<<","<<pressureDivY<<","<<stressDivY
//                 <<endl;
//                 exit(1);
//             }
//             else
//                 PP->V[n] += ((2.0/3.0)*dv2 + (2.0/3.0)*dv3)*p->dt;
//             if(dw2!=dw2||dw3!=dw3)
//             {
//                 cerr<<"Particle velocity component w resulted in NaN.\n"
//                 <<dw2<<","<<dw3<<"|"<<Dp<<","<<netBuoyZ<<","<<pressureDivZ<<","<<stressDivZ
//                 <<endl;
//                 exit(1);
//             }
//             else
//                 PP->W[n] += ((2.0/3.0)*dw2 + (2.0/3.0)*dw3)*p->dt;

//             // if(print)
//             // {
//             //     cout<<((2.0/3.0)*du2 + (2.0/3.0)*du3)<<","<<((2.0/3.0)*dv2 + (2.0/3.0)*dv3)<<","<<((2.0/3.0)*dw2 + (2.0/3.0)*dw3)<<"\n"
//             //     <<Dp*du<<","<<netBuoyX<<","<<-pressureDivX/p->S22<<","<<-stressDivX/((1-thetas)*p->S22)<<"\n"
//             //     <<Dp*dv<<","<<netBuoyY<<","<<-pressureDivY/p->S22<<","<<-stressDivY/((1-thetas)*p->S22)<<"\n"
//             //     <<Dp*dw<<","<<netBuoyZ<<","<<-pressureDivZ/p->S22<<","<<-stressDivZ/((1-thetas)*p->S22)<<"\n"
//             //     <<stressDivZ<<","<<thetas<<"\n"
//             //     <<p->mpirank<<"|"<<i<<","<<j<<","<<k
//             //     <<endl;
//             //     print=false;
//             // }
            
//             // Pos update

//             // Solid forcing
//             double solid_old = p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n]);
//             double solid_new = p->ccipol4_b(a->solid,PP->X[n]+PP->U[n]*p->dt,PP->Y[n]+PP->V[n]*p->dt,PP->Z[n]+PP->W[n]*p->dt);
//             if(solid_new<=0)
//             {
//                 double solid_x = p->ccipol4_b(a->solid,PP->X[n]+PP->U[n]*p->dt,PP->Y[n],PP->Z[n]);
//                 double solid_y = p->ccipol4_b(a->solid,PP->X[n],PP->Y[n]+PP->V[n]*p->dt,PP->Z[n]);
//                 double solid_z = p->ccipol4_b(a->solid,PP->X[n],PP->Y[n],PP->Z[n]+PP->W[n]*p->dt);
//                 if(solid_x<=0)
//                 {
//                     double dx = (solid_old)/(solid_old-solid_x)*PP->U[n]*p->dt+(PP->U[n]>=0?-1:1)*PP->d50/2.0;
//                     PP->X[n] += dx;
//                     PP->U[n] = 0;
//                 }
//                 if(solid_y<=0)
//                 {
//                     double dy = (solid_old)/(solid_old-solid_y)*PP->V[n]*p->dt+(PP->V[n]>=0?-1:1)*PP->d50/2.0;
//                     PP->Y[n] += dy;
//                     PP->W[n] = 0;
//                 }
//                 if(solid_z<=0)
//                 {
//                     double dz = (solid_old)/(solid_old-solid_z)*PP->W[n]*p->dt+(PP->W[n]>=0?-1:1)*PP->d50/2.0;
//                     PP->Z[n] += dz;
//                     PP->W[n] = 0;
//                 }
//             }

//             PP->X[n] += PP->U[n]*p->dt;
//             PP->Y[n] += PP->V[n]*p->dt;
//             PP->Z[n] += PP->W[n]*p->dt;

//             // Sum update
//             cellSum[IJK]-=PP->PackingFactor[n];
//             if(cellSum[IJK]<0)
//             {
//                 clog<<"cellSum is below zero in cell ("<<p->XN[IP]<<"-"<<p->XN[IP1]<<","<<p->YN[JP]<<"-"<<p->YN[JP1]<<","<<p->ZN[KP]<<"-"<<p->ZN[KP1]<<") of partition "<<p->mpirank<<" for particle "<<n<<" at ("<<PP->X[n]<<","<<PP->Y[n]<<","<<PP->Z[n]<<") "<<solid_new<<"."<<endl;
//                 cellSum[IJK]=0;
//             }
//             particleStressTensorUpdateIJK(p,a,PP);
//             i=p->posc_i(PP->X[n]);
//             j=p->posc_j(PP->Y[n]);
//             k=p->posc_k(PP->Z[n]);
//             cellSum[IJK]+=PP->PackingFactor[n];
//             particleStressTensorUpdateIJK(p,a,PP);
//         }
// }