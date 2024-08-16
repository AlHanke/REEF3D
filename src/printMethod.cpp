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

#include "printMethod.h"

#include "lexer.h"

printMethod::printMethod(lexer *p)
{
    switch (p->P10)
    {
        case 0:
            outputFormat = new vtk3D();
            break;
        case 1: default:
            outputFormat = new vtu3D();
            break;
        case 2:
            outputFormat = new vtr3D();
            break;
        case 3:
            outputFormat = new vts3D();
            break;
    }

    if(p->mpirank==0)
    {
        outputFormat->folder("CFD");
    }
}
printMethod::~printMethod()
{
    delete outputFormat;
}