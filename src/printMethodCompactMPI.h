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

#ifndef PRINT_METHOD_COMPACT_MPI_H_
#define PRINT_METHOD_COMPACT_MPI_H_

#include "printMethod.h"

#include "field5.h"
#include "nodefill.h"

#include <mpi.h>

#include <sstream>
#include <vector>

class printMethodCompactMPI : public printMethod
{
    public:
        printMethodCompactMPI(lexer*);
        ~printMethodCompactMPI()=default;
        void setup(lexer*,fdm*,ghostcell*);
        int print(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
    private:
        void offsetCMPIPoints(lexer*,int*,int*,int);
    private:
        size_t headerSize=0;
        int endIndex=0;
        std::vector<MPI_Offset> offsetCMPI;
        std::vector<int> offsetCMPIitr;
        std::stringstream header;
        std::stringstream footer;
        int kbeginPoint,kendPoint;
        int jbeginPoint,jendPoint;
        int ibeginPoint,iendPoint;
        int compactMPIPOffset[300];
        int cellNum;
        int pointNum;

        field5 eta;
        nodefill node;

        int m;
};

#endif