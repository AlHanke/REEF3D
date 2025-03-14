/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"force.h"
#include<stdio.h>

void force::pvtp(int M10)
{   
    sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%06i.pvtp",ID);

    ofstream result;
    result.open(name);

    result<<"<?xml version=\"1.0\"?>\n";
    result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    result<<"<PPolyData GhostLevel=\"0\">\n";
    
    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>\n";
    result<<"</PPointData>\n";

    result<<"<PPoints>\n";
    result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    result<<"</PPoints>\n";

    for(n=0; n<M10; ++n)
    {
        snprintf(pname,sizeof(pname),"REEF3D-SOLID-%06i-%06i.vtp",ID,n+1);
        result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    result<<"</PPolyData>\n";
    result<<"</VTKFile>"<<endl;

    result.close();
}
