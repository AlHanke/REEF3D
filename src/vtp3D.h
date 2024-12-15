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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef VTP3D_H_
#define VTP3D_H_

#include "vtk3D.h"
#include "increment.h"

class lexer;

class vtp3D : public vtk3D , virtual public increment
{
    protected:
        void beginning(lexer*, std::stringstream&, int, int, int, int, int);
        void beginningParallel(lexer*, std::ofstream&);

        void points(std::stringstream&, const int*, int&);
        void pointsParallel(std::ofstream&);
        void verts(std::stringstream&, const int*, int&);
        void polys(std::stringstream&, const int*, int&);

        void ending(std::stringstream&);
        void endingParallel(std::ofstream&);
        void footer(std::vector<char>&, int&);
};

#endif