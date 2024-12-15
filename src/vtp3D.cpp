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

#include "vtp3D.h"

#include "lexer.h"

void vtp3D::beginning(lexer *p, std::stringstream &result, int numPoints, int numVerts, int numLines, int numStrips, int numPolys)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PolyData\" ";
    vtkVersion(result);
    result<<"<PolyData>\n";
    if(p->P16==1)
        timeValue(result,p->simtime);
    result<<"<Piece NumberOfPoints=\""<<numPoints<<"\" NumberOfVerts=\""<<numVerts<<"\" NumberOfLines=\""<<numLines<<"\" NumberOfStrips=\""<<numStrips<<"\" NumberOfPolys=\""<<numPolys<<"\">\n";
}

void vtp3D::beginning(lexer *p, std::ofstream &result)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PolyData\" ";
    vtkVersion(result);
    result<<"<PolyData>\n";
    if(p->P16==1)
        timeValue(result,p->simtime);
}

void vtp3D::beginningParallel(lexer *p, std::ofstream &result)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PPolyData\" ";
    vtkVersion(result);
    result<<"<PPolyData GhostLevel=\"0\">\n";
    if(p->P16==1)
        timeValue(result,p->simtime);
}