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

#include"topo_vtp.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sstream>
#include<cstdio>
#include<cstring>

void topo_vtp::print(lexer* p, fdm* a, ghostcell *pgc, sediment *psed)
{
    if(p->mpirank==0)
        pvtp(p,a,pgc,psed);
    
    name_iter(p,a,pgc);

    //---------------------------------------------
    n=0;
    offset[n]=0;
    ++n;
    
    //Velocity
    offset[n]=offset[n-1] + 4*p->pointnum2D*3+ 4;
    ++n;
    // Elevation
    offset[n]=offset[n-1] + 4*p->pointnum2D + 4;
    ++n;
    
    // sediment bedlaod
    if(p->P76==1)
        psed->offset_ParaView_2D_bedload(p,pgc,offset,n);

    // sediment parameters 1
    if(p->P77==1)
        psed->offset_ParaView_2D_parameter1(p,pgc,offset,n);

    // sediment parameters 2
    if(p->P78==1)
        psed->offset_ParaView_2D_parameter2(p,pgc,offset,n);

    // bed shear stress
    if(p->P79>=1)
        psed->offset_ParaView_2D_bedshear(p,pgc,offset,n);
    
    //End Data

    // Points
    offset[n]=offset[n-1] + 4*p->pointnum2D*3 + 4;
    ++n;
    
    // Polys
    // Connectivities
    offset[n]=offset[n-1] + 4*polygon_sum*3 + 4;
    ++n;
    // Offsets
    offset[n]=offset[n-1] + 4*polygon_sum + 4;
    ++n;
    //---------------------------------------------

    std::stringstream result;

    beginning(p,result,p->pointnum2D,0,0,0,polygon_sum);

    n=0;
    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    
    if(p->P76==1)
        psed->name_ParaView_bedload(p,pgc,result,offset,n);
    
    if(p->P77==1)
        psed->name_ParaView_parameter1(p,pgc,result,offset,n);

    if(p->P78==1)
        psed->name_ParaView_parameter2(p,pgc,result,offset,n);

    if(p->P79>=1)
        psed->name_ParaView_bedshear(p,pgc,result,offset,n);
    result<<"</PointData>\n";

    points(result,offset,n);

    polys(result,offset,n);

    ending(result);
    
    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);

    //----------------------------------------------------------------------------
    
    //  Velocities
    iin=4*(p->pointnum2D)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol1a(a->P));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(p->sl_ipol2a(a->Q));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=0.0;
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    //  Elevation
    iin=4*p->pointnum2D;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
    {
        ffn=float(p->sl_ipol4(a->bed));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    //  sediment bedload
    if(p->P76==1)
        psed->print_2D_bedload(p,pgc,buffer,m);
    
    //  sediment parameter 1
    if(p->P77==1)
        psed->print_2D_parameter1(p,pgc,buffer,m);

    //  sediment parameter 2
    if(p->P78==1)
        psed->print_2D_parameter2(p,pgc,buffer,m);

    //  bed shear stress
    if(p->P79>=1)
        psed->print_2D_bedshear(p,pgc,buffer,m);

    //  Points XYZ
    iin=4*(p->pointnum2D)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    TPSLICELOOP
    {
        ffn=p->XN[IP1];
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=p->XN[IP1];
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=p->XN[IP1];
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    // Polys

    //  Connectivity
    iin=4*(polygon_sum)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    SLICEBASELOOP
    {
        // Triangle 1
        iin=int(a->nodeval2D(i-1,j-1))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);

        iin=int(a->nodeval2D(i,j-1))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);

        iin=int(a->nodeval2D(i,j))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);


        // Triangle 2
        iin=int(a->nodeval2D(i-1,j-1))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);

        iin=int(a->nodeval2D(i,j))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);

        iin=int(a->nodeval2D(i-1,j))-1;
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
    }


    // Offset of Connectivity
    iin=4*(polygon_sum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<polygon_sum;++n)
    {
    iin=(n+1)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    }

    footer(buffer,m);

    // Open File
    FILE* file = std::fopen(name, "w");
    std::fwrite(buffer.data(), buffer.size(), 1, file);
    std::fclose(file);
}







