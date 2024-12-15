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

#include"particle_pls.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void particle_pls::pvtu_pos(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	snprintf(name,sizeof(name),"./REEF3D_PLS/XPLS-POS-%08i.pvtu",num);


	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<PUnstructuredGrid GhostLevel=\"0\">\n";

	result<<"<PPointData>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"radius\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"correction\"/>\n";
	result<<"</PPointData>\n";

	result<<"<PPoints>\n";
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
	result<<"</PPoints>\n";

	for(n=0; n<p->M10; ++n)
	{
    piecename_pos(p,a,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	result<<"</PUnstructuredGrid>\n";
	result<<"</VTKFile>\n";

	result.close();
}


void particle_pls::pvtu_neg(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	
	snprintf(name,sizeof(name),"./REEF3D_PLS/XPLS-NEG-%08i.pvtu",num);
	
	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	result<<"<PUnstructuredGrid GhostLevel=\"0\">\n";

	result<<"<PPointData>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"radius\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"correction\"/>\n";
	result<<"</PPointData>\n";

	result<<"<PPoints>\n";
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
	result<<"</PPoints>\n";

	for(n=0; n<p->M10; ++n)
	{
    piecename_neg(p,a,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	result<<"</PUnstructuredGrid>\n";
	result<<"</VTKFile>\n";

	result.close();
}

void particle_pls::piecename_pos(lexer* p, fdm* a, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;


    snprintf(pname,sizeof(pname),"XPLS-POS-%08i-%08i.vtu",num,n+1);

}

void particle_pls::piecename_neg(lexer* p, fdm* a, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	snprintf(pname,sizeof(pname),"XPLS-NEG-%08i-%08i.vtu",num,n+1);
}

void particle_pls::header_pos(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

    
    snprintf(name,sizeof(name),"./REEF3D_PLS/XPLS-POS-%08i-%06i.vtp",num,p->mpirank+1);

}

void particle_pls::header_neg(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;


    snprintf(name,sizeof(name),"./REEF3D_PLS/XPLS-NEG-%08i-%06i.vtp",num,p->mpirank+1);

}
