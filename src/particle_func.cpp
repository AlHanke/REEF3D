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

#include"particle_func.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"tracers_obj.h"
#include"particles_obj.h"
#include"boundarycheck.h"

#define PARTICLELOOP for(size_t n=0;n<PP->loopindex;n++)

particle_func::particle_func(lexer* p) : kinVis(p->W1/p->W2), drho(p->W1/p->S22)
{
    p->Darray(stressTensor,p->imax*p->jmax*p->kmax);
}
particle_func::~particle_func()
{
}

/// @brief Applies advection to positions of particles in @param PP
/// @param p partition object
/// @param a fdm object contains flow field
/// @param PP tracers_obj contains tracer information
/// @param minflag PP.Flag[n] needs to be bigger than minflag for PP[n] to be affected by avection
void particle_func::advect(lexer* p, fdm* a, tracers_obj* PP, int minflag, double source_u, double source_v, double source_w)
{
    double coord1, coord2, coord3, u1, u2, v1, v2, w1, w2;
    PARTICLELOOP
        if(PP->Flag[n]>minflag)
        {
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            PP->cellSum[IJK]--;

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
            PP->cellSum[IJK]++;
        }
}

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
            PP->cellSum[IJK]--;

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
            PP->cellSum[IJK]++;
        }
}

/// @brief Particle transport equation following (Tavouktsoglou 2021)
/// @param p 
/// @param a 
/// @param PP 
/// @param minflag Particles need to have a larger flag
/// Uses Runge Kutta 3 to calculate the change in particle velocity 
/// The particle velocity is then used to transport the particle
void particle_func::transport(lexer* p, fdm* a, particles_obj* PP, int minflag)
{
    double RKu,RKv,RKw;
    double u,v,w;
    double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
    /// @brief Difference between flowfield and particle velocity
    double du, dv, dw;
    double Dp, thetas;
    double pressureDiv=0, stressDiv=0;

    PARTICLELOOP
        if(PP->Flag[n]>minflag)
        {
            // Prep
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);

            thetas=0;//(PI*pow(PP->d50,3.0)*PP->cellSum[IJK]/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]));
            u=p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n]);
            v=p->ccipol1(a->v,PP->X[n],PP->Y[n],PP->Z[n]);
            w=p->ccipol1(a->w,PP->X[n],PP->Y[n],PP->Z[n]);

            // RK3 step 1
            du=u-PP->U[n];
            dv=v-PP->V[n];
            dw=w-PP->W[n];

            Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

            du1=Dp*du+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dv1=Dp*dv+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dw1=Dp*dw+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));

            RKu=PP->U[n]+du1*p->dt;
            RKv=PP->V[n]+dv1*p->dt;
            RKw=PP->W[n]+dw1*p->dt;
            
            // RK step 2
            du=u-RKu;
            dv=v-RKv;
            dw=w-RKw;

            Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

            du2=Dp*du+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dv2=Dp*dv+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dw2=Dp*dw+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));

            du2=0.25*du2+0.25*du1;
            dv2=0.25*dv2+0.25*dv1;
            dw2=0.25*dw2+0.25*dw1;

            RKu=PP->U[n]+du2*p->dt;
            RKv=PP->V[n]+dv2*p->dt;
            RKw=PP->W[n]+dw2*p->dt;
            
            // RK step 3
            du=u-RKu;
            dv=v-RKv;
            dw=w-RKw;

            Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

            du3=Dp*du+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dv3=Dp*dv+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
            dw3=Dp*dw+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));

            PP->U[n] += ((2.0/3.0)*du2 + (2.0/3.0)*du3)*p->dt;
            PP->V[n] += ((2.0/3.0)*dv2 + (2.0/3.0)*dv3)*p->dt;
            PP->W[n] += ((2.0/3.0)*dw2 + (2.0/3.0)*dw3)*p->dt;

            // Pos update
            PP->X[n] += PP->U[n]*p->dt;
            PP->Y[n] += PP->V[n]*p->dt;
            PP->Z[n] += PP->W[n]*p->dt;

            // Sum update
            PP->cellSum[IJK]--;
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            PP->cellSum[IJK]++;
        }
}

/// @brief Removes all tracers with have left the partions boundaries
/// @param p partition object
/// @param PP tracers_obj contains tracer information
/// @return Number of removed tracers
int particle_func::remove(lexer* p,tracers_obj* PP)
{
    bool inBounds=false;
    int removed=0;
    int i,j,k;
    boundarycheck bounderies;

    PARTICLELOOP
        if(PP->Flag[n]>0)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);

            inBounds=bounderies.minboundcheck(p,i,j,k,1);
            if (inBounds)
                inBounds=bounderies.maxboundcheck(p,i,j,k,1);

			// remove out of bounds particles
            if(!inBounds)
            {
                PP->erase(n);
                PP->cellSum[IJK]--;
                removed++;
            }
        }
    
    return removed;
}

/// @brief Transfer function
/** Is responsible to transfer tracers to one of the surrounding partitions*/
/// @param p partition object
/// @param pgc ghostcell object
/// @param PP tracers_obj contains tracer information
/// @param maxcount maximum number of tracers which could be transfered
/// @return number of send of tracers
int particle_func::transfer(lexer* p, ghostcell* pgc, tracers_obj* PP, int maxcount)
{
    int xchange=0;

    tracers_obj seedling1(maxcount),seedling2(maxcount),seedling3(maxcount),seedling4(maxcount),seedling5(maxcount),seedling6(maxcount);
    tracers_obj Send[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};
    tracers_obj Recv[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};

    int i,j,k;

    PARTICLELOOP
        if(PP->Flag[n]==1)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);


            if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
            {
                switch (p->flag5[IJK])
                {
                    case -1:
                    {
                        Send[0].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -2:
                    {
                        Send[1].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -3:
                    {
                        Send[2].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -4:
                    {
                        Send[3].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -5:
                    {
                        Send[4].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }

                    case -6:
                    {
                        Send[5].add(PP->X[n],PP->Y[n],PP->Z[n],1);
                        break;
                    }
                }
                PP->erase(n);
                PP->cellSum[IJK]--;
                ++xchange;
            }
        }

    pgc->para_tracersobj(p,Send,Recv);


    size_t sum=PP->size;
    for(int n=0;n<6;n++)
        sum += Recv[n].size;
    if(sum>PP->capacity)
        PP->reserve(sum);

    for(int n=0;n<6;n++)
        PP->add_obj(&Recv[n]);
    
    return xchange;
}


/// @brief Transfer function
/// Is responsible to transfer particles to one of the surrounding partitions
/// @param p partition object
/// @param pgc ghostcell object
/// @param PP particles_obj contains particle information
/// @param maxcount maximum number of particles which could be transfered
/// @return number of send of particles
int particle_func::transfer(lexer* p, ghostcell* pgc, particles_obj* PP, int maxcount)
{
    int xchange=0;

    particles_obj seedling1(maxcount,PP->d50,PP->density,true),seedling2(maxcount,PP->d50,PP->density,true),seedling3(maxcount,PP->d50,PP->density,true),seedling4(maxcount,PP->d50,PP->density,true),seedling5(maxcount,PP->d50,PP->density,true),seedling6(maxcount,PP->d50,PP->density,true);
    particles_obj Send[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};
    particles_obj Recv[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};

    int i,j,k;

    PARTICLELOOP
        if(PP->Flag[n]==1)
        {
            i = p->posc_i(PP->X[n]);
            j = p->posc_j(PP->Y[n]);
            k = p->posc_k(PP->Z[n]);


            if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
            {
                switch (p->flag5[IJK])
                {
                    case -1:
                    {
                        Send[0].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -2:
                    {
                        Send[1].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -3:
                    {
                        Send[2].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -4:
                    {
                        Send[3].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -5:
                    {
                        Send[4].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }

                    case -6:
                    {
                        Send[5].add(PP->X[n],PP->Y[n],PP->Z[n],1,PP->U[n],PP->V[n],PP->W[n],PP->PackingFactor[n]);
                        break;
                    }
                }
                PP->erase(n);
                PP->cellSum[IJK]--;
                ++xchange;
            }
        }

    pgc->para_tracersobj(p,Send,Recv);


    size_t sum=PP->size;
    for(int n=0;n<6;n++)
        sum += Recv[n].size;
    if(sum>PP->capacity)
        PP->reserve(sum);

    for(int n=0;n<6;n++)
        PP->add_obj(&Recv[n]);

    return xchange;
}

/// @brief Particle Reynolds number
/// Calculates particle reynolds number for particle[index]
/// @param p
/// @param a
/// @param PP 
/// @param index
/// @return Local particle Reynolds number
double particle_func::reynolds(lexer* p,fdm* a, particles_obj* PP, int index)
{
    const double u=p->ccipol1(a->u,PP->X[index],PP->Y[index],PP->Z[index]);
    const double v=p->ccipol2(a->v,PP->X[index],PP->Y[index],PP->Z[index]);
    const double w=p->ccipol3(a->w,PP->X[index],PP->Y[index],PP->Z[index]);

    const double mean_vel=sqrt(u*u+v*v+w*w);

    const double Re=mean_vel*PP->d50/p->W2;
    // Change to particle diameter once implemented

    return Re;
}

/// @brief Settling velocity
/// Calculates settling velocity using drag coefficent
/// g is assumed to act only in negative z direction
/// @param p 
/// @param PP 
/// @return 
double particle_func::settling_vel(lexer* p,fdm* a, particles_obj* PP, int index)
{
    return sqrt(4.0/3.0*(PP->density/p->W1-1.0)*fabs(p->W22)*PP->d50/drag_coefficient(p,a,PP,index));
}

/// @brief Drag coefficent Cd
/// Calculates drag coefficent from particle Reynolds number based on emperical 
/// @return 
double particle_func::drag_coefficient(lexer* p,fdm* a, particles_obj* PP, int index)
{
    const double Re=reynolds(p,a,PP,index);
    if(Re<0.5)
        return 24.0/Re;
    if(Re<1000.0)
        return 18.5*pow(Re,-0.6);
    if(Re<200000.0)
        return 0.44;
    return NAN;


    /// Andrews and O’Rourke (1996)
    const double Rep=0;
    const double theta_s=0;// vol sol/vol cell
    const double theta_f=1-theta_s;
    return 24.0*(pow(theta_f,-2.65)+pow(Rep,2.0/3.0)*pow(theta_f,-1.78)/6.0)/Rep;
}

/// @brief 
/// @param p 
/// @param PP 
void particle_func::make_stationary(lexer* p, fdm* a, tracers_obj* PP, int minflag)
{
    int i,j;
    PARTICLELOOP
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])<0)
            PP->Flag[n]=minflag;
}

/// @brief 
/// @param p 
/// @param PP 
void particle_func::make_stationary(lexer* p, fdm* a, particles_obj* PP, int minflag)
{
    int i,j;
    PARTICLELOOP
        if(PP->Flag[n]>0)
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])<0)
        {
            PP->Flag[n]=minflag;
            if(p->count!=0)
            {
                i=p->posc_i(PP->X[n]);
                j=p->posc_j(PP->Y[n]);
                p->flag_topo_changed[IJ]=1;
                p->topo_change[IJ]+=volume(PP,n);
            }
            if(PP->entries>PP->tracers_obj::entries)
            {
                PP->U[n]=0;
                PP->V[n]=0;
                PP->W[n]=0;
            }
        }
}

/// @brief 
/// @param PP 
/// @param index 
/// @return 
double particle_func::volume(particles_obj* PP, int index)
{
    return PI*pow(PP->d50,3.0)/6.0;
}

void particle_func::cleanup(lexer* p, fdm* a, particles_obj* PP, int max)
{
    int* numPartijk;
    int i,j,k;

    p->Iarray(numPartijk,p->knox*p->knoy*p->knoz);
    PARTICLELOOP
        if(PP->Flag[n]==0)
        {
            i=p->posc_i(PP->X[n]);
            j=p->posc_j(PP->Y[n]);
            k=p->posc_k(PP->Z[n]);
            if(a->topo(i,j,k)<p->DZN[KP])
                if(numPartijk[IJK]<max)
                    numPartijk[IJK]++;
                else
                    PP->erase(n);
        }
    p->del_Iarray(numPartijk,p->knox*p->knoy*p->knoz);
}

/// @brief Topo volume in cell div. by particle volume
/// Uses i,j&k from increment to pass cell identifier
/// @param d50 Sauter diameter of particles
/// @return Ceil of number of particles in cell IJK
int particle_func::maxParticlesPerCell(lexer* p, fdm* a, double d50)
{   
    return ceil(p->DXN[IP]*p->DYN[JP]*(a->topo(i,j,k)<p->DZN[KP]?a->topo(i,j,k):p->DZN[KP])*(1-p->S24)*6.0/(PI*pow(d50,3.0)));
}

int particle_func::maxParticlesPerXY(lexer* p, fdm* a, double d50)
{
    return ceil(p->DXN[IP]*p->DYN[JP]*(1-p->S24)*6.0/(PI*pow(d50,2.0)));
}

void particle_func::particlesPerCell(lexer* p, particles_obj* PP)
{
    PARTICLELOOP
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        k=p->posc_k(PP->Z[n]);
        PP->cellSum[IJK]++;
    }
}

void particle_func::particleStressTensor(lexer* p, particles_obj* PP)
{
    double theta;
    int i,j,k;

    PLAINLOOP
    {
        theta = PI*pow(PP->d50,3.0)*PP->cellSum[IJK]/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
        stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1-theta));
    }
}

void particle_func::particleStressTensorUpdateIJK(lexer* p, particles_obj* PP)
{
    double theta;
    int i,j,k;

    for (int n=-2; n<3; n++)
        for (int m=-2; m<3; m++)
            for (int l=-2; l<3; l++)
            {
                i=increment::i+n;
                j=increment::j+m;
                k=increment::k+l;

                theta = PI*pow(PP->d50,3.0)*PP->cellSum[IJK]/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]);
                stressTensor[IJK]=Ps*pow(theta,beta)/max(theta_crit-theta,epsilon*(1-theta));
            }
}

double particle_func::drag_model(lexer* p, double d, double du, double dv, double dw, double thetas) const
{
    const double thetaf = 1.0-thetas;

    const double dU=sqrt(du*du+dv*dv+dw*dw);

    const double Rep=dU*d*kinVis;

    const double Cd=24.0*(pow(thetaf,-2.65)+pow(Rep,2.0/3.0)*pow(thetaf,-1.78)/6.0)/Rep;
    const double Dp=Cd*3.0*drho*dU/d/4.0;

    return Dp;
}

void particle_func::make_moving(lexer* p, fdm* a, particles_obj* PP)
{
    double RKu,RKv,RKw;
    double u,v,w;
    double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
    /// @brief Difference between flowfield and particle velocity
    double du, dv, dw;
    double Dp, thetas;
    double pressureDiv=0, stressDiv=0;

    PARTICLELOOP
    if(PP->Flag[n]==0)
    {
        i=p->posc_i(PP->X[n]);
        j=p->posc_j(PP->Y[n]);
        k=p->posc_k(PP->Z[n]);

        thetas=0;//(PI*pow(PP->d50,3.0)*PP->cellSum[IJK]/(6.0*p->DXN[IP]*p->DYN[JP]*p->DYN[KP]));
        u=p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n]);
        v=p->ccipol1(a->v,PP->X[n],PP->Y[n],PP->Z[n]);
        w=p->ccipol1(a->w,PP->X[n],PP->Y[n],PP->Z[n]);

        // RK3 step 1
        du=u-PP->U[n];
        dv=v-PP->V[n];
        dw=w-PP->W[n];

        Dp=drag_model(p,PP->d50,du,dv,dw,thetas);

        du1=Dp*du+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
        dv1=Dp*dv+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));
        dw1=Dp*dw+(1.0-drho)*p->W20-(pressureDiv/p->S22+stressDiv/((1-thetas)*p->S22));

        if (fabs(du1)>0||fabs(dv1)>0||dw1>0)
        PP->Flag[n]=1;
    }
}