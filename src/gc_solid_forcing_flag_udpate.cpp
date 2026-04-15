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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing_flag_update(lexer *p, fdm *a)
{
    // Update DF
    set_DF(p, a);
    
    // 1
    ULOOP
    {
        int df=p->DF(i,j,k);

        if(df>0 && p->DF(i+1,j,k)<0)
            df=-1;

        p->DF1[IJK]=df;
    }

    // 2
    VLOOP
    {
        int df=p->DF(i,j,k);

        if(df>0 && p->DF(i,j+1,k)<0)
            df=-1;

        p->DF2[IJK]=df;
    }

    // 3
    WLOOP
    {
        int df=p->DF(i,j,k);

        if(df>0 && p->DF(i,j,k+1)<0)
            df=-1;

        p->DF3[IJK]=df;
    }

    #if USE_AMREX
    ArrayWrapper3D::FillBoundaryBatch(p, p->m_df123, 0, 3);
    #else
    gcparaxintV(p, p->DF1, 1);
    gcparaxintV(p, p->DF2, 1);
    gcparaxintV(p, p->DF3, 1);
    #endif
}
