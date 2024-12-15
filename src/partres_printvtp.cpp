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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include<sstream>
#include<cstdio>
#include<cstring>

void partres::print_particles(lexer* p, sediment_fdm *s)
{
    if((p->count%p->Q181==0 || p->count==0) && (p->Q180==1 && p->Q181>0 && p->Q182<0.0))
    {
    print_vtp(p,s);
    ++printcount;
    }
    
    if((p->simtime>p->partprinttime || p->count==0) && (p->Q180==1 && p->Q181<0 && p->Q182>0.0))
    {
    print_vtp(p,s);
    p->partprinttime+=p->Q182;
    ++printcount;
    }
}

void partres::print_vtp(lexer* p, sediment_fdm *s)
{
    int numpt=0;

    for(n=0;n<P.index;++n)
    if(P.Flag[n]>0)
    numpt++;

    //cout<<"PSed-"<<p->mpirank<<"| printed: "<<numpt<<" not printed: "<<P.size-numpt<<" | capcacity: "<<P.capacity<<endl;

    int count;
    int n=0;
    int offset[100];
    int iin;
    float ffn;
    
    if(p->mpirank==0)
    pvtp(p);

    header_pos(p);

    offset[n]=0;
    ++n;
    
    offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //flag
    ++n;
    if(p->P23==1)
    {
        offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //Test
        ++n;
    }
    offset[n]=offset[n-1]+sizeof(float)*(numpt)*3+sizeof(int); //velocity
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //radius
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(numpt)*3+sizeof(int); //fluid velocity
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //bedChange
    ++n;

    offset[n]=offset[n-1]+sizeof(float)*(numpt)*3+sizeof(int); //xyz
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //connectivitey
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*(numpt)+sizeof(int); //offset connectivity
    ++n;

    //---------------------------------------------
    stringstream result;

    beginning(p,result);
    result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfVerts=\""<<numpt<<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    
    result<<"<PointData>\n";
    n=0;
    result<<"<DataArray type=\"Int32\" Name=\"Flag\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P23==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"Test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"fluid velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"bedChange\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";
    

    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";
    

    result<<"<Verts>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Verts>\n";

    result<<"</Piece>\n";
    result<<"</PolyData>\n";

    result<<"<AppendedData encoding=\"raw\">\n_";
    m=result.str().length();
    buffer.resize(m+offset[n]+27);
    std::memcpy(&buffer[0],result.str().data(),m);
    //----------------------------------------------------------------------------
    
    // flag
    iin=sizeof(int)*(numpt);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        iin=int(p->mpirank);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
    }
    
    // Test
    if(p->P23==1)
    {
        iin=sizeof(float)*(numpt);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.Test[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }
    }

    // velocities
    iin=sizeof(float)*(numpt)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        ffn=float(P.U[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.V[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.W[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    // radius
    iin=sizeof(float)*(numpt);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        ffn=float(P.d50/2);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    // fluid velocities
    iin=sizeof(float)*(numpt)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        ffn=float(P.Uf[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.Vf[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.Wf[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    // bedChange
    iin=sizeof(float)*(numpt);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        //ffn=float(p->ccslipol4(s->bedch,P.X[n],P.Y[n]));
        ffn=float(p->ccslipol4(s->bedzh,P.X[n],P.Y[n])-p->ccslipol4(s->bedzh0,P.X[n],P.Y[n]));
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }

    //  XYZ
    iin=sizeof(float)*(numpt)*3;
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        ffn=float(P.X[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.Y[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);

        ffn=float(P.Z[n]);
        std::memcpy(&buffer[m],&ffn,sizeof(float));
        m+=sizeof(float);
    }
    
    //  Connectivity
    iin=sizeof(int)*(numpt);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    count=0;
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        iin=int(count);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        ++count;
    }

    //  Offset of Connectivity
    iin=sizeof(int)*(numpt);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    count=1;
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
    {
        iin=int(count);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        ++count;
    }

    stringstream footer;
    footer<<"\n</AppendedData>\n</VTKFile>"<<flush;
    std::memcpy(&buffer[m],footer.str().data(),footer.str().size());

    // Open File
    FILE* file = fopen(name, "w");
    fwrite(buffer.data(), buffer.size(), 1, file);
    fclose(file);
}
