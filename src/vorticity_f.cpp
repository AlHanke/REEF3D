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

#include"vorticity_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sstream>
#include<vector>
#include<cstring>

vorticity_f::vorticity_f(lexer *p, fdm *a) : strain(p,a), omega1(p), omega2(p), omega3(p)
{
}

void vorticity_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, int &m)
{
    double wx,wy,wz;

// xy plane

    LOOP
    {
     omega1(i,j,k) = qij(p,a,2,3);
     omega2(i,j,k) = qij(p,a,1,3);
     omega3(i,j,k) = qij(p,a,2,1);
    }

    pgc->start4(p,omega1,1);
	pgc->start4(p,omega2,1);
	pgc->start4(p,omega3,1);

    // --
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
	{
    ffn=float(p->ipol4(omega1));

	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
	
	// --
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
	{
    ffn=float(p->ipol4(omega2));

	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
	
	// --
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
	{
    ffn=float(p->ipol4(omega3));

	std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
	}
}

void vorticity_f::name_ParaView_parallel(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"vorticity x\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"vorticity y\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"vorticity z\"/>\n";
}

void vorticity_f::name_ParaView(lexer *p, fdm *a, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"vorticity x\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"vorticity y\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"vorticity z\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void vorticity_f::offset_ParaView(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}




