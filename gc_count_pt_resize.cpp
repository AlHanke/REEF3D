/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::column_pt_resize(lexer* p, fdm* a)
{
    int size=0;
    int safety=1000;
    
    size = column_pt1_count(p,a);
    
    if(size>p->C1_size)
    {
    cout<<p->mpirank<<" CPT1 Resize: "<<p->C1_size<<" "<<size<<endl;
    size+=safety;
    a->C1.resize(p,p->C1_size,size);
    
    p->C1_size=size;
    }

    
    
    size = column_pt2_count(p,a);
    
    if(size>p->C2_size)
    {
    cout<<p->mpirank<<" CPT2 Resize: "<<p->C2_size<<" "<<size<<endl;
    size+=safety;
    a->C2.resize(p,p->C2_size,size);
    
    p->C2_size=size;
    }
    
    
    size = column_pt3_count(p,a);
    
    if(size>p->C3_size)
    {
    cout<<p->mpirank<<" CPT3 Resize: "<<p->C3_size<<" "<<size<<endl;
    size+=safety;
    a->C3.resize(p,p->C3_size,size);
    
    p->C3_size=size;
    }
    
    
    
    size = column_pt4_count(p,a);
    
    if(size>p->C4_size)
    {
    cout<<p->mpirank<<" CPT4 Resize: "<<p->C4_size<<" "<<size<<endl;
    size+=safety;
    a->C4.resize(p,p->C4_size,size);
    
    p->C4_size=size;
    }
    
    
    size = column_pt4a_count(p,a);
    
    if(size>p->C4a_size)
    {
    cout<<p->mpirank<<" CPT4a Resize: "<<p->C4a_size<<" "<<size<<endl;
    size+=safety;
    a->C4a.resize(p,p->C4a_size,size);
    
    p->C4a_size=size;
    }
    
    
    size = column_pt6_count(p,a);
    
    if(size>p->C6_size)
    {
    cout<<p->mpirank<<" CPT6 Resize: "<<p->C6_size<<" "<<size<<endl;
    size+=safety;
    a->C6.resize(p,p->C6_size,size);
    
    p->C6_size=size;
    }
   
}