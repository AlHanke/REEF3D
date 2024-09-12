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

#include"ghostcell.h"

#include "lexer.h"

#include "tracers_obj2.h"
#include "particles_obj2.h"


void ghostcell::setup_MPI_particle_dataTypes(lexer* p)
{
    {
        const int number_of_blocks = 4;
        int blocklengths[number_of_blocks] = {1, 1, 1, 1};
        MPI_Datatype types[number_of_blocks] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
        MPI_Aint offsets[number_of_blocks];
        offsets[0] = offsetof(tracers_obj2::tracer, x);
        offsets[1] = offsetof(tracers_obj2::tracer, y);
        offsets[2] = offsetof(tracers_obj2::tracer, z);
        offsets[3] = offsetof(tracers_obj2::tracer, flag);
        
        MPI_Type_create_struct(number_of_blocks, blocklengths, offsets, types, &MPI_TRACER);
        MPI_Type_commit(&MPI_TRACER);
    }

    {
        const int number_of_blocks = 11;
        int blocklengths[number_of_blocks] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        MPI_Datatype types[number_of_blocks] = {MPI_TRACER, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        MPI_Aint offsets[number_of_blocks];
        offsets[0] = 0;
        offsets[1] = offsetof(particles_obj2::particle, u);
        offsets[2] = offsetof(particles_obj2::particle, v);
        offsets[3] = offsetof(particles_obj2::particle, w);
        offsets[4] = offsetof(particles_obj2::particle, parcelFactor);
        offsets[5] = offsetof(particles_obj2::particle, uf);
        offsets[6] = offsetof(particles_obj2::particle, vf);
        offsets[7] = offsetof(particles_obj2::particle, wf);
        offsets[8] = offsetof(particles_obj2::particle, shear_eff);
        offsets[9] = offsetof(particles_obj2::particle, shear_crit);
        offsets[10] = offsetof(particles_obj2::particle, drag);
        
        MPI_Type_create_struct(number_of_blocks, blocklengths, offsets, types, &MPI_PARTICLE);
        MPI_Type_commit(&MPI_PARTICLE);
    }
    
}

void ghostcell::free_MPI_particle_dataTypes()
{
    MPI_Type_free(&MPI_TRACER);
    MPI_Type_free(&MPI_PARTICLE);
}