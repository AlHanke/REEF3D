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
#include"particles_obj.h"

class lexer;
class fdm;
class ghostcell;

class tracers_obj;
//class particles_obj;

/// Particle function class
/** A class containing all basic function to manipulate the position of tracers_objs. */
class particle_func: private increment
{
public:

protected:
    particle_func(lexer*,int=10,double=0.001,double=2700);
    virtual ~particle_func();
    
    // para
    int remove(lexer*,tracers_obj*);
    int remove(lexer*,particles_obj*);
    int transfer(lexer*,ghostcell*,tracers_obj*,int);
    int transfer(lexer*,ghostcell*,particles_obj*,int);

    // mov
    void advect(lexer*,fdm*,tracers_obj*,int=0,double=0,double=0,double=0);
    void advect(lexer*,fdm*,particles_obj*,int=0,double=0,double=0,double=0);
    void transport(lexer*,fdm*,particles_obj*,int=0);
    void make_stationary(lexer*,fdm*,tracers_obj*,int=0);
    void make_stationary(lexer*,fdm*,particles_obj*,int=0);
    void make_moving(lexer*,fdm*,particles_obj*);

    // util
    double reynolds(lexer*,fdm*,particles_obj*,int);
    double settling_vel(lexer*,fdm*,particles_obj*,int);
    double drag_coefficient(lexer*,fdm*,particles_obj*,int);
    double volume(particles_obj*,int);
    double maxParticlesPerCell(lexer*,fdm*,double,bool=true);
    int maxParticlesPerXY(lexer*,fdm*,double);
    void particlesPerCell(lexer*,ghostcell*,particles_obj*);
    void particleStressTensor(lexer*,fdm*,ghostcell*,particles_obj*);
    void particleStressTensorUpdateIJK(lexer*,fdm*,particles_obj*);
    void updateParticleStressTensor(lexer*,fdm*,particles_obj*,int,int,int);
    double theta_s(lexer*,fdm*,particles_obj*,int,int,int);
    double drag_model(lexer*,double,double,double,double,double) const;
    void debug(lexer*,fdm*,ghostcell*,particles_obj*);
    void fixPos(lexer*,fdm*,particles_obj*);

    // memory management
    void cleanup(lexer*,fdm*,particles_obj*,int);
private:
    const double kinVis;
    const double drho;
    double* stressTensor;
    double* cellSum;
    const double Ps; // in pressure unit
    const double beta; // 2<=beta<=5
    const double epsilon;
    const double theta_crit; // 0.6-0.65
    particles_obj seedling1,seedling2,seedling3,seedling4,seedling5,seedling6;
};

#endif