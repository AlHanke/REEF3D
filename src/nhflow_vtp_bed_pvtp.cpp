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

#include"nhflow_vtp_bed.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment.h"
#include<stdio.h>

void nhflow_vtp_bed::pvtp(lexer *p, fdm_nhf *d, ghostcell* pgc, sediment *psed)
{    
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
    
    std::snprintf(name,sizeof(name),"./REEF3D_NHFLOW_VTP_BED/REEF3D-NHFLOW-BED-%08i.pvtp",num);


    ofstream result;
    result.open(name);

    beginningParallel(p,result);
    
    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"depth\"/>\n";
    
    if(p->P76==1)
    psed->name_ParaView_parallel_bedload(p,pgc,result);
    
    if(p->P77==1)
    psed->name_ParaView_parallel_parameter1(p,pgc,result);

    if(p->P78==1)
    psed->name_ParaView_parallel_parameter2(p,pgc,result);

    if(p->P79>=1)
    psed->name_ParaView_parallel_bedshear(p,pgc,result);
    
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";
    result<<"</PPointData>\n";

    pointsParallel(result);

    for(n=0; n<p->M10; ++n)
    {
    piecename(p,d,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    endingParallel(result);

    result.close();

}

void nhflow_vtp_bed::piecename(lexer *p, fdm_nhf *d, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

    std::snprintf(pname,sizeof(pname),"REEF3D-NHFLOW-BED-%08i-%06i.vtp",num,n+1);

}
