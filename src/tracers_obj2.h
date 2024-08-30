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

#ifndef TRACERSOBJ2_H_
#define TRACERSOBJ2_H_

#include <stdio.h>
#include <cstdint>
#include <vector>

class tracers_obj2
{
public:
    tracers_obj2(size_t _capacity=10, size_t _size=0);
    virtual ~tracers_obj2()=default;

    virtual void erase(size_t _index);
    void erase_all();
    size_t add(double x, double y, double z, int flag);
    void add_obj(tracers_obj2* obj);
    size_t add_entry(tracers_obj2* obj, size_t _index);
    void print(size_t _index);
public:
    struct tracer
    {
        double x;
        double y;
        double z;
        int flag;
    };
    std::vector<tracer> tracers;
};

#endif
