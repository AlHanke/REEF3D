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

#ifndef PARTICLESOBJ2_H_
#define PARTICLESOBJ2_H_

#include "tracers_obj2.h"

class lexer;

class particles_obj2 : public tracers_obj2
{
public:

    struct particle : tracers_obj2::tracer
    {
        double u;
        double v;
        double w;
        double parcelFactor=1;
        double xrk1=0;
        double yrk1=0;
        double zrk1=0;
        double urk1=0;
        double vrk1=0;
        double wrk1=0;
        double uf=0;
        double vf=0;
        double wf=0;
        double shear_eff=0;
        double shear_crit=0;
        double drag=0;
    };

public:
    particles_obj2(size_t=10, double=0.001, double=2650.0);
    virtual ~particles_obj2()=default;

    void erase(size_t);
    void erase_all();

    size_t add_entry(particle);
    void add_obj(particles_obj2*);

    void print(particle);

public:
    // --- particle data ---
    // ---    general    ---

    /// @brief d50 of particle set
    double d50;
    /// @brief Average density of particle set
    const double density;

    std::vector<particle> particles;
};

#endif
