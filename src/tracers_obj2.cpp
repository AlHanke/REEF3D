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

#include "tracers_obj2.h"

#include <iostream>

/// @brief Tracer style particles\n
/// Contains and manages massless particles used as tracers in the flow field
/// @param capacity Desired initial capacity
/// @param size Desired number of tracers at default position (0,0,0|INT32_MIN)
/// @param scale_factor Sets ::scale_factor for ::reserve
tracers_obj2::tracers_obj2(size_t _capacity, size_t _size)
{	
    tracers.reserve(_capacity);
    for(size_t n=0;n<_size;n++)
    {
        tracer temp;
        temp.x=0;
        temp.y=0;
        temp.z=0;
        temp.flag=INT32_MIN;
        tracers.push_back(temp);
    }
}

/// @brief Removes \p _index
void tracers_obj2::erase(size_t _index)
{
    tracers.erase(tracers.begin()+_index);
}

/// @brief Clears all data
void tracers_obj2::erase_all()
{
    tracers.clear();
}

/// @brief Addes new particle with prescribed position and state
/// @param x Position in x-dir
/// @param y Position in y-dir
/// @param z Position in z-dir
/// @param flag State - stationary, moving, etc.
/// @return Index of added particle
size_t tracers_obj2::add(double x, double y, double z, int flag)
{
    tracers.push_back({x,y,z,flag});
    return tracers.size()-1;
}

/// @brief Addes new element based on \p _index of \p obj
/// @return Index of newly added element
size_t tracers_obj2::add_entry(tracers_obj2* obj, size_t _index)
{
    tracers.push_back(obj->tracers[_index]);
    return tracers.size()-1;
}

/// @brief Adds contens of object of type tracers_obj
/// @param obj Inputed data
void tracers_obj2::add_obj(tracers_obj2* obj)
{
    tracers.insert(tracers.end(),obj->tracers.begin(),obj->tracers.end());
}

/// @brief Prints position of element to cout
/// @param _index Element to print
void tracers_obj2::print(size_t _index)
{
    std::cout<<"Tracer_obj["<<_index<<"]=("<<tracers[_index].x<<","<<tracers[_index].y<<","<<tracers[_index].y<<")"<<std::endl;
}
