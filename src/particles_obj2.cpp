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
#include <iostream>

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
particles_obj2::particles_obj2(size_t _capacity, double _d50, double _density):
                d50(_d50), density(_density)
                
{	
    particles.reserve(_capacity);
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

/// \copydoc tracers_obj::add_entry
size_t particles_obj2::add_entry(particle _particle)
{
    particles.push_back(_particle);
    return particles.size()-1;
}

/// \copydoc tracers::add_obj
void particles_obj2::add_obj(particles_obj2* obj)
{
    if(obj->density!=density)
        std::cerr<<"particles_obj::add_obj - density mismatch"<<std::endl;
    if(obj->d50!=d50)
        std::cerr<<"particles_obj::add_obj - d50 mismatch"<<std::endl;
    // for (auto particle : obj->particles)
    //     std::cout<<particle.x<<" "<<particle.y<<" "<<particle.z<<std::endl;
    
    particles.insert(particles.end(),obj->particles.begin(),obj->particles.end());
}

void particles_obj2::print(particle _particle)
{
    std::cout<<"Sending\n"<<_particle.x<<" "<<_particle.y<<" "<<_particle.z<<"\n"
    <<_particle.xrk1<<" "<<_particle.yrk1<<" "<<_particle.zrk1<<"\n"
    <<_particle.u<<" "<<_particle.v<<" "<<_particle.w<<"\n"
    <<_particle.urk1<<" "<<_particle.vrk1<<" "<<_particle.wrk1<<"\n"<<std::endl;
}