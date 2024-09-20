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

#ifndef PRINT_METHOD_COMPACT_H_
#define PRINT_METHOD_COMPACT_H_

#include "printMethod.h"

class printMethodCompact : public printMethod
{
    public:
        printMethodCompact(lexer*);
        ~printMethodCompact();
        void setup(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
        int print(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
    private:
        double *XN=nullptr;
        double *YN=nullptr;
        double *ZN=nullptr;

        double *  uvelGlobal=nullptr;
        double ** uvel=nullptr;
        double *  vvelGlobal=nullptr;
        double ** vvel=nullptr;
        double *  wvelGlobal=nullptr;
        double ** wvel=nullptr;

        double *  pressGlobal=nullptr;
        double ** press=nullptr;
        
        double *  eddyvGlobal=nullptr;
        double ** eddyv=nullptr;
        
        double *  phiGlobal=nullptr;
        double ** phi=nullptr;

        double *  rhoGlobal=nullptr;
        double ** rho=nullptr;

        double *  viscGlobal=nullptr;
        double ** visc=nullptr;

        double *  VOFGlobal=nullptr;
        double ** VOF=nullptr;

        double *  topoGlobal=nullptr;
        double ** topo=nullptr;

        double *  testGlobal=nullptr;
        double ** test=nullptr;

        double *  solidGlobal=nullptr;
        double ** solid=nullptr;

        double *  fbGlobal=nullptr;
        double ** fb=nullptr;

        double *  walldGlobal=nullptr;
        double ** walld=nullptr;

        
        int *  flagGlobal=nullptr;
        int ** flag=nullptr;
        
        int *  flag5Global=nullptr;
        int ** flag5=nullptr;
        
        
        int *globalSendCounts=nullptr;
        int *displs=nullptr;
        
        int localSendCount=0;
        int cellNum=0;
        int pointNum=0;
};

#endif