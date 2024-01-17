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

class lexer;
class fdm;
class ghostcell;

class tracers_obj;

/// Particle function class
/** A class containing all basic function to manipulate the position of tracers_objs. */
class particle_func
{
public:
    particle_func();
    virtual ~particle_func(); 
    void advect(lexer*,fdm*,tracers_obj*,int);
    int remove(lexer*,tracers_obj*);
    int transfer(lexer*,ghostcell*,tracers_obj*,int);
private:

};

#endif