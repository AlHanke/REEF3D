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

#include "particles_obj2.h"
#include"lexer.h"
#include <cstdint>
#include <cstring>
#include<iostream>

/*
Dangers when using:
size_t overflow when adding something to an object at capacity
*/

/// @brief Container for particles with density, d50 and x,y,z velocities
/// @param _capacity Desired initial capacity
/// @param _d50 \copydoc particles_obj::d50
/// @param _density \copydoc particles_obj::density
/// @param individuals Whether particles are allowed to have individual data besides position
/// @param _size Desired number of partices at default position (0,0,0|INT32_MIN) with individual data (0,0,0|0)
/// @param _scale_factor Sets ::scale_factor for ::reserve
particles_obj2::particles_obj2(size_t _capacity, double _d50, double _density, size_t _size):
                d50(_d50), density(_density)
                
{	
    particles.reserve(_capacity);
    for(size_t n=0;n<_size;n++)
    {
        particle temp;
        temp.x=0;
        temp.y=0;
        temp.z=0;
        temp.flag=INT32_MIN;
        temp.u=0;
        temp.v=0;
        temp.w=0;
        temp.packingFactor=1;
        temp.uf=0;
        temp.vf=0;
        temp.wf=0;
        temp.shear_eff=0;
        temp.shear_crit=0;
        temp.drag=0;
        particles.push_back(temp);
    }
}

/// \copydoc tracers_obj::erase
void particles_obj2::erase(size_t _index)
{
    particles.erase(particles.begin()+_index);
}

/// \copydoc tracers_obj::erase_all
void particles_obj2::erase_all()
{
    particles.clear();
}

/// @brief Addes new particle with prescribed position and state
/// @param x Position in x-dir
/// @param y Position in y-dir
/// @param z Position in z-dir
/// @param flag State - stationary, moving, etc.
/// @param u Velocity in x-dir
/// @param v Velocity in y-dir
/// @param w Velocity in z-dir
/// @param packingFactor Number of real particles represented by the element
/// @return Index of added particle
size_t particles_obj2::add(const double x, const double y, const double z, const int flag, const double u, const double v, const double w, const double packingFactor, const double uF, const double vF, const double wF, const double shearEff, const double shearCrit, const double _drag)
{
    // particles.push_back({x,y,z,flag,u,v,w,packingFactor,uF,vF,wF,shearEff,shearCrit,_drag});
    // particles.push_back({{x, y, z, flag}, u, v, w, packingFactor, uF, vF, wF, shearEff, shearCrit, _drag});
    particle temp;
    temp.x=x;
    temp.y=y;
    temp.z=z;
    temp.flag=flag;
    temp.u=u;
    temp.v=v;
    temp.w=w;
    temp.packingFactor=packingFactor;
    temp.uf=uF;
    temp.vf=vF;
    temp.wf=wF;
    temp.shear_eff=shearEff;
    temp.shear_crit=shearCrit;
    temp.drag=_drag;
    particles.push_back(temp);
    return particles.size()-1;
}

/// \copydoc tracers_obj::add_entry
size_t particles_obj2::add_entry(particles_obj2* obj, size_t _index)
{
    particles.push_back(obj->particles[_index]);
    return particles.size()-1;
}

/// \copydoc tracers::add_obj
void particles_obj2::add_obj(particles_obj2* obj)
{
    if(obj->density!=density)
        std::cerr<<"particles_obj::add_obj - density mismatch"<<std::endl;
    if(obj->d50!=d50)
        std::cerr<<"particles_obj::add_obj - d50 mismatch"<<std::endl;
    
    particles.insert(particles.end(),obj->particles.begin(),obj->particles.end());
}
