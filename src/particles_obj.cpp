/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"particles_obj.h"
#include"lexer.h"
#include <cstdint>
#include <cstring>
#include<iostream>

/*
Dangers when using:
size_t overflow when adding something to an object at capacity
*/


particles_obj::particles_obj(size_t capacity, double d50, double density, bool individuals, size_t size, double scale_factor):
                d50(d50), density(density), scale_factor(scale_factor), tracers_obj(capacity,size,scale_factor),
                flag_inactive(0), flag_bed(1), flag_bed_load(2), flag_suspended_load(3),
                entries(tracers_obj::entries+(individuals?4:0)) // update when adding more data
{	
    if(capacity>0)
    {
        if(size>capacity)
            capacity=size;

        if(individuals)
        {
            U = new double[capacity];
            V = new double[capacity];
            W = new double[capacity];
            
            PackingFactor = new double[capacity];
        }
    }
}

particles_obj::~particles_obj()
{
}

/// \copydoc tracers_obj::debug
void particles_obj::debug()
{
    std::cout<<"particle_obj::debug"<<std::endl;
    // insert code for debugging here //
    
}

/// \copydoc tracers_obj::erase
void particles_obj::erase(size_t index)
{
    tracers_obj::erase(index);

    if(entries>tracers_obj::entries)
    {
        U[index]=0;
        V[index]=0;
        W[index]=0;

        PackingFactor[index]=0;
    }
}

/// \copydoc tracers_obj::erase_all
void particles_obj::erase_all()
{
    tracers_obj::erase_all();

    if(entries>tracers_obj::entries)
    {
        delete[] U;
        delete[] V;
        delete[] W;
        U = new double[capacity];
        V = new double[capacity];
        W = new double[capacity];

        delete[] PackingFactor;
        PackingFactor = new double[capacity];
    }
}

/// @brief
/// @param x 
/// @param y 
/// @param z 
/// @param flag 
/// @param u vel
/// @param v 
/// @param w 
/// @param packingFactor Number of particles in parcel
/// @return 
size_t particles_obj::add(double x, double y, double z, int flag, double u, double v, double w, double packingFactor)
{
    std::cout<<size<<"|"<<capacity<<"|"<<loopindex<<"|"<<Empty[empty_itr]<<std::endl;
    size_t index=tracers_obj::add(x,y,z,flag);
    // std::cout<<size<<"|"<<capacity<<"|"<<loopindex<<"|"<<index<<std::endl;
    if(entries>tracers_obj::entries)
        add_data(index,u,v,w,packingFactor);
    return index;
}

/// \copydoc tracers_obj::reserve
size_t particles_obj::reserve(size_t capacity_desired)
{
    if(0==capacity_desired)
        capacity_desired=ceil(scale_factor*capacity);
    if (capacity_desired>capacity)
    {
        if (capacity_desired>SIZE_MAX)
            std::__throw_length_error("particles_obj - max capacity reached");

        tracers_obj::reserve(capacity_desired);
    }
    return this->capacity;
}

/// \copydoc tracers_obj::fill
void particles_obj::fill(size_t index, bool do_empty, int flag)
{
    if(entries>tracers_obj::entries)
        for(size_t n=size; n<index;++n)
        {
            U[n]=0;
            V[n]=0;
            W[n]=0;

            PackingFactor[n]=1;
        }
    tracers_obj::fill(index,do_empty,flag);
}

/// \copydoc tracers_obj::check_state
bool particles_obj::check_state(bool first)
{
    if(capacity-size!=empty_itr||size>capacity)
    {
        if(first)
        {
            fix_state();
            return check_state(false);
        }
        else
            return false;
    }
    else
        return true;
}

/// \copydoc tracers_obj::fix_state
void particles_obj::fix_state() // ToDo - update
{
    if(this->size>this->capacity)
    {
        size_t real_size=capacity;
        size_t old_size=size;
        reserve(ceil(this->scale_factor*size));
        this->size=real_size;
        fill(old_size,false);

        this->empty_itr=0;
        size_t temp_Empty[capacity];
        for(size_t n=0;n<real_size;n++)
            if(Empty[n]!=NULL)
                temp_Empty[empty_itr++]=n;
        delete [] Empty;
        this->Empty=temp_Empty;
        fill_empty();
    }
    else
    {
        this->empty_itr=0;
        size_t temp_Empty[capacity];
        delete [] Empty;
        this->Empty=temp_Empty;
        for(size_t n=capacity-1; n>0;--n)
            if(X[n]==NULL&&Y[n]==NULL&&Z[n]==NULL&&Flag[n]==-1)
                Empty[empty_itr++]=n;
    }
}

/// \copydoc tracers_obj::optimize
void particles_obj::optimize() // ToDo - update
{
    // Could be optimized to look ahead if a large section is empty to only do one move operation
    size_t loopchange=0;
    for(int n=empty_itr; n>=0;--n)
        if(Empty[n]<loopindex)
        {
            memorymove(Empty[n],Empty[n]+1,(loopindex-Empty[n]+1));

            loopchange++;
        }
    loopindex -= loopchange;
}

/// \copydoc tracers::memorymove
void particles_obj::memorymove(size_t des, size_t src, size_t len)
{
    tracers_obj::memorymove(des,src,len);

    if(entries>tracers_obj::entries)
    {
        std::memmove(&U[des],&U[src],sizeof(double)*len);
        std::memmove(&V[des],&V[src],sizeof(double)*len);
        std::memmove(&W[des],&W[src],sizeof(double)*len);

        std::memmove(&PackingFactor[des],&PackingFactor[src],sizeof(double)*len);
    }
}

/// @brief 
/// @param obj 
void particles_obj::add_obj(particles_obj* obj)
{
    if(obj->size>0)
    {
        if(size+obj->size>capacity)
            reserve(size+obj->size);
        
        if(obj->entries>obj->tracers_obj::entries && this->entries>this->tracers_obj::entries)
            for(size_t n=0;n<obj->loopindex;n++)
                add(obj->X[n],obj->Y[n],obj->Z[n],obj->Flag[n],obj->U[n],obj->V[n],obj->W[n],obj->PackingFactor[n]);
        else
            tracers_obj::add_obj(obj);
    }
}

/// \copydoc tracers_obj::add_obj
void particles_obj::add_obj(tracers_obj* obj)
{
    tracers_obj::add_obj(obj);
}

/// @brief Additional data input
/// @param index 
/// @param u 
/// @param v 
/// @param w 
/// @param packingFactor 
void particles_obj::add_data(size_t index, double u, double v, double w, double packingFactor)
{
    std::cout<<size<<"|"<<capacity<<"|"<<loopindex<<"|"<<index<<std::endl;
    U[index]=u;
    V[index]=v;
    W[index]=w;
    // PackingFactor[index]=packingFactor;
}