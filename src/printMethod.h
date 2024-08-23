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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef PRINT_METHOD_H_
#define PRINT_METHOD_H_

#include "increment.h"
#include "vtks.h"

#include <fstream>

class lexer;
class fdm;
class ghostcell;

class print_averaging;
class turbulence;
class heat;
class multiphase;
class vorticity;
class data;
class concentration;
class sediment;

class printMethod : protected virtual increment
{
    public:
        printMethod(lexer *);
        virtual ~printMethod();
        virtual void setup(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*) = 0;
        virtual int print(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*) = 0;
    protected:
        void calcVTKOffsets(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*,const int,const int);
    
    protected:
        vtk3D *outputFormat=nullptr;
        char name[200];
        int vtkOffsets[300];
        int n=0;
        int iin=0;
        float ffn=0;
};

#endif