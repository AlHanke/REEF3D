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

#ifndef PRINT_METHOD_SEPARATED_H_
#define PRINT_METHOD_SEPARATED_H_

#include "printMethod.h"

#include "field5.h"
#include "nodefill.h"

class printMethodSeparated : public printMethod
{
    public:
        printMethodSeparated(lexer*);
        ~printMethodSeparated()=default;
        void setup(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
        int print(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
    private:
        void parallelData(lexer*,fdm*,ghostcell*,print_averaging*,turbulence*,heat*,multiphase*,vorticity*,data*,concentration*,sediment*);
    
    private:
        field5 eta;
        nodefill node;
        int filecount=0;
};

#endif