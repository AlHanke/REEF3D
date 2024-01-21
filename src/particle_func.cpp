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

particle_func::particle_func()
{
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
    for(size_t n=0;n<PP->loopindex;n++)
        if(PP->Flag[n]>minflag)
        {
            u1=p->dt*p->ccipol1(a->u,PP->X[n],PP->Y[n],PP->Z[n]);
            coord1=PP->X[n]+u1;
            
            v1=p->dt*p->ccipol2(a->v,PP->X[n],PP->Y[n],PP->Z[n]);
            coord2=PP->Y[n]+v1;
            
            w1=p->dt*p->ccipol3(a->w,PP->X[n],PP->Y[n],PP->Z[n]);
            coord3=PP->Z[n]+w1;
            
            
            u2=0.25*u1 + 0.25*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);
            coord1=PP->X[n]+u2;
            
            v2=0.25*v1 + 0.25*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
            coord2=PP->Y[n]+v2;
            
            w2=0.25*w1 + 0.25*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);
            coord3=PP->Z[n]+w2;
            
            
            PP->X[n] += (2.0/3.0)*u2 + (2.0/3.0)*p->dt*p->ccipol1(a->u,coord1,coord2,coord3);

            PP->Y[n] += (2.0/3.0)*v2 + (2.0/3.0)*p->dt*p->ccipol2(a->v,coord1,coord2,coord3);
            
            PP->Z[n] += (2.0/3.0)*w2 + (2.0/3.0)*p->dt*p->ccipol3(a->w,coord1,coord2,coord3);
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

    for(size_t n=0;n<PP->loopindex;n++)
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

    for(size_t n=0;n<PP->loopindex;n++)
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
    {
        PP->add_obj(&Recv[n]);
        // Send[n].erase_all();
        // Recv[n].erase_all();
    }
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
    return sqrt(4.0/3.0*(PP->density/p->W2-1.0)*fabs(p->W22)*PP->d50/drag_coefficient(p,a,PP,index));
}

/// @brief Drag coefficent Cd
/// Calculates drag coefficent from particle Reynolds number based on emperical 
/// @return 
double particle_func::drag_coefficient(lexer* p,fdm* a, particles_obj* PP, int index)
{
    const double Re=reynolds(p,a,PP,index);
    if(Re<0.5)
        return 24.0/Re;
    if(Re<1000)
        return 18.5*pow(Re,-0.6);
    if(Re<200000)
        return 0.44;
    return NAN;
}

/// @brief 
/// @param p 
/// @param PP 
void particle_func::make_stationary(lexer* p, fdm* a, tracers_obj* PP)
{
    int i,j;
    for(size_t n=0;n<PP->loopindex;n++)
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])<0)
            PP->Flag[n]=0;
}

/// @brief 
/// @param p 
/// @param PP 
void particle_func::make_stationary(lexer* p, fdm* a, particles_obj* PP)
{
    int i,j;
    for(size_t n=0;n<PP->loopindex;n++)
        if(PP->Flag[n]>0)
        if (p->ccipol4_b(a->topo,PP->X[n],PP->Y[n],PP->Z[n])<0)
        {
            PP->Flag[n]=0;
            if(p->count!=0)
            {
                i=p->posc_i(PP->X[n]);
                j=p->posc_j(PP->Y[n]);
                p->flag_topo_changed[IJ]=1;
                p->topo_change[IJ]+=volume(PP,n);
                cout<<"Topo increased by particles in cell("<<i<<","<<j<<")"<<IJ<<endl;
            }
        }
}

/// @brief 
/// @param PP 
/// @param index 
/// @return 
double particle_func::volume(particles_obj* PP, int index)
{
    return PI*PP->d50*PP->d50*PP->d50/6;
}