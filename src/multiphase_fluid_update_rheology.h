/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef MULTIPHASE_FLUID_UPDATE_RHEOLOGY_H_
#define MULTIPHASE_FLUID_UPDATE_RHEOLOGY_H_

#include"multiphase_fluid_update.h"
#include"increment.h"

class rheology;

using namespace std;

class multiphase_fluid_update_rheology : public multiphase_fluid_update, increment
{
public:
    multiphase_fluid_update_rheology(lexer*);
    virtual ~multiphase_fluid_update_rheology();

    void start(lexer*, fdm*, ghostcell*,field&,field&) override;

private:
    const int gcval_ro,gcval_visc;
    const double ro1;
    double visc1;
    const double ro2,visc2;
    const double ro3,visc3;
    rheology *prheo;
};

#endif
