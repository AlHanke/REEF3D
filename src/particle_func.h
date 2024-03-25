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
    
    // Parallelization
    int remove(lexer*,tracers_obj*);
    int transfer(lexer*,ghostcell*,tracers_obj*,int);
    int transfer(lexer*,ghostcell*,particles_obj*,int);

    // Movement
    void advect(lexer*,fdm*,tracers_obj*,int=0,double=0,double=0,double=0);
    void advect(lexer*,fdm*,particles_obj*,int=0,double=0,double=0,double=0);
    void transport(lexer*,fdm*,particles_obj*,int=0);
    void make_stationary(lexer*,fdm*,tracers_obj*,int=0);
    void make_stationary(lexer*,fdm*,particles_obj*,int=0);
    void make_moving(lexer*,fdm*,particles_obj*);

    // Utility
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
    void debug(lexer*,fdm*,ghostcell*,tracers_obj*);
    void fixPos(lexer*,fdm*,tracers_obj*);

    // memory management
    void cleanup(lexer*,fdm*,tracers_obj*,int);
protected:
    /// @brief Inter-particle stresses per cell
    double* stressTensor;
    /// @brief Number of particles in a cell
    double* cellSum;
    /// @brief Volume change per column
    double* topoVolumeChange;

private:
    /// @brief Kinetic viscory of fluid\n Initialized using `lexer::W1` and `lexer::W2`
    const double kinVis;
    /// @brief Ratio of fluid and solid densities\n Initialized using `lexer::W1` and `lexer::S22`
    const double drho;
    /// @brief Constant for stress trensor calculation\n Initialized using `lexer::Q14` and given in Pascal
    const double Ps;
    /// @brief Constant for stress trensor calculation\n Initialized using `lexer::Q15` and should be in range of \f$2\leq\beta\leq5\f$
    const double beta;
    /// @brief Dampener for stress trensor calculation\n Used to dampen out sudden acceleration resulting from high packing densities, initialized using `lexer::Q16`
    const double epsilon;
    /// @brief Maximum solid volume fraction of a fully packed bed\n Usually between 60% and 65%
    const double theta_crit;
};

#endif