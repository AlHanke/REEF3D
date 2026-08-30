/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fieldint4.h"

void ghostcell::gcsldf_update(lexer *p)
{
    count=0;
    
    // eta
    k=p->knoz-1;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        ++count;
        
        if(p->DF[Ip1JK]<0)
        ++count;
        
        if(p->DF[IJm1K]<0)
        ++count;
        
        if(p->DF[IJp1K]<0)
        ++count;
    }
    
    
    if(p->gcsldfeta4_count!=count)
    {
        p->gcsldfeta4_count=count;

        p->gcsldfeta4.resize(p->gcsldfeta4_count);
    }
    
    
    
    count=0;
    k=p->knoz-1;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        {
        p->gcsldfeta4[count][0]=i;
        p->gcsldfeta4[count][1]=j;
        p->gcsldfeta4[count][2]=1;
        ++count;
        }
        
        if(p->DF[Ip1JK]<0)
        {
        p->gcsldfeta4[count][0]=i;
        p->gcsldfeta4[count][1]=j;
        p->gcsldfeta4[count][2]=4;
        ++count;
        }
        
        if(p->DF[IJm1K]<0)
        {
        p->gcsldfeta4[count][0]=i;
        p->gcsldfeta4[count][1]=j;
        p->gcsldfeta4[count][2]=3;
        ++count;
        }
        
        if(p->DF[IJp1K]<0)
        {
        p->gcsldfeta4[count][0]=i;
        p->gcsldfeta4[count][1]=j;
        p->gcsldfeta4[count][2]=2;
        ++count;
        }
    }
    
    // -------------------------
    // bed
    k=0;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        ++count;
        
        if(p->DF[Ip1JK]<0)
        ++count;
        
        if(p->DF[IJm1K]<0)
        ++count;
        
        if(p->DF[IJp1K]<0)
        ++count;
    }
    
    
    if(p->gcsldfbed4_count!=count)
    {
        p->gcsldfbed4_count=count;
        p->gcsldfbed4.resize(p->gcsldfbed4_count);
    }
    
    
    // assign gcsldfbed entries
    
    count=0;
    
    k=0;
    SLICELOOP4
    if(p->DF[IJK]==1)
    {
        if(p->DF[Im1JK]<0)
        {
        p->gcsldfbed4[count][0]=i;
        p->gcsldfbed4[count][1]=j;
        p->gcsldfbed4[count][2]=1;
        ++count;
        }
        
        if(p->DF[Ip1JK]<0)
        {
        p->gcsldfbed4[count][0]=i;
        p->gcsldfbed4[count][1]=j;
        p->gcsldfbed4[count][2]=4;
        ++count;
        }
        
        if(p->DF[IJm1K]<0)
        {
        p->gcsldfbed4[count][0]=i;
        p->gcsldfbed4[count][1]=j;
        p->gcsldfbed4[count][2]=3;
        ++count;
        }
        
        if(p->DF[IJp1K]<0)
        {
        p->gcsldfbed4[count][0]=i;
        p->gcsldfbed4[count][1]=j;
        p->gcsldfbed4[count][2]=2;
        ++count;
        }
    }
    
    
    
    
}