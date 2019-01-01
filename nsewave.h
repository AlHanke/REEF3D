/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/


class convection;
class pressure;
class solver;
class fdm;
class lexer;
class ghostcell;
class fluid_update;
class heat;
class concentration;
class ioflow;
class slice;
class momentum;
class diffusion;
class poisson;
class turbulence;

using namespace std;

#ifndef NSEWAVE_H_
#define NSEWAVE_H_

class nsewave 
{
public:    
    virtual void start(lexer*, fdm*, ghostcell*, momentum*, diffusion*, turbulence*, convection*, 
                        pressure*, poisson*, solver*, solver*, ioflow*)=0;
    
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*)=0;

        

};

#endif
