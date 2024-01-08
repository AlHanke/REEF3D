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

#ifndef PARTICLESOBJ_H_
#define PARTICLESOBJ_H_

#include <stdio.h>
#include <math.h>

/*
Philosophy: performance, memory usage, ease of use

All data should be stored as at max a 1D array
Internal data access should be via iterator
No external access to iterator

Out of bounds safe
Thread safe
*/

class particles_obj
{
public:
    particles_obj(size_t, double=0.001, double=2650.0, double=0.4, size_t=0, double=1.25);
    ~particles_obj();

    void erase(size_t);
    void reserve(size_t);
    void add(double,double,double,int);
    bool check_state(bool=true);

private:
    void fill(size_t,bool=true);
    void fill_empty();
    void fix_state();
    void debug();
    void clear();

public:
    // state data
    size_t size;
    size_t capacity;

    //particle data
    //general
    double d50;
    const double density;
    const double porosity;

    //individual
    double* X;
    double* Y;
    double* Z;

    int* Flag;

private:
    size_t empty_itr;
    size_t* Empty;
    const double scale_factor;
};

#endif
