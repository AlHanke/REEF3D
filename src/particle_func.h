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

#ifndef PARTICLEFUNC_H_
#define PARTICLEFUNC_H_

#include"increment.h"

class lexer;
class fdm;
class ghostcell;

class tracers_obj;
class particles_obj;

/// Particle function class
/** A class containing all basic function to manipulate the position of tracers_objs. */
class particle_func: private increment
{
public:

protected:
    particle_func();
    virtual ~particle_func();
    
    // para
    int remove(lexer*,tracers_obj*);
    int transfer(lexer*,ghostcell*,tracers_obj*,int);
    int transfer(lexer*,ghostcell*,particles_obj*,int);

    // mov
    void advect(lexer*,fdm*,tracers_obj*,int=0,double=0,double=0,double=0);
    void advect(lexer*,fdm*,particles_obj*,int=0,double=0,double=0,double=0);
    void transport(lexer*,fdm*,particles_obj*,int=0);
    void make_stationary(lexer*,fdm*,tracers_obj*);
    void make_stationary(lexer*,fdm*,particles_obj*);

    // util
    double reynolds(lexer*,fdm*,particles_obj*,int);
    double settling_vel(lexer*,fdm*,particles_obj*,int);
    double drag_coefficient(lexer*,fdm*,particles_obj*,int);
    double volume(particles_obj*,int);
    int maxParticlesPerCell(lexer*,fdm*,double);
    int maxParticlesPerXY(lexer*,fdm*,double);
    void particlesPerCell(lexer*,particles_obj*);

    // memory management
    void cleanup(lexer*,fdm*,particles_obj*,int);
private:
    const double kinVis;
    const double drho;
    double* stressTensor;
    double Ps=1; // in pressure unit
    double beta = 3.5; // 2<=beta<=5
    double epsilon = 10e-7;
    double theta_crit = 0.6; // 0.6-0.65
};

#endif