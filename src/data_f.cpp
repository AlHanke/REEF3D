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

#include"data_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cstring>

data_f::data_f(lexer* p, fdm *a, ghostcell* pgc) : data(p)
{
}

data_f::~data_f()
{
}

void data_f::start(lexer* p, fdm* a, ghostcell* pgc)
{
	cout<<"DATA "<<p->P150<<endl;
	
	if(p->P151==1)
    LOOP
    data(i,j,k)=p->data[(i-p->imin)*p->jmax + (j-p->jmin)];
	
	if(p->P151==2)
	LOOP
	data(i,j,k)=-p->data[(i-p->imin)*p->jmax + (j-p->jmin)]+p->pos_z();
	
	if(p->P152==1)
	pgc->start4(p,data,101);
	
	if(p->P152==2)
	pgc->start4(p,data,102);
	
	if(p->P152==3)
	pgc->start4(p,data,103);
	
	if(p->P152==4)
	pgc->start4(p,data,1);
}

void data_f::print_3D(lexer* p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, int &m)
{
    iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

    TPLOOP
    {
    ffn=float(p->ipol4(data));
    std::memcpy(&buffer[m],&ffn,sizeof(float));
    m+=sizeof(float);
    }
}

void data_f::name_ParaView_parallel(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"data\"/>\n";
}

void data_f::name_ParaView(lexer *p, fdm *a, ghostcell *pgc, stringstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"data\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void data_f::offset_ParaView(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}


