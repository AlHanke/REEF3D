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

void ghostcell::flagbase(lexer *p)
{
    int bc=0;

    MALOOP
        p->flag5[IJK]=0;

    GC4LOOP
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        if(p->gcb4[n][4]==1)
        bc=1;
        else if(p->gcb4[n][4]==2)
        bc=2;
        else if(p->gcb4[n][4]==3)
        bc=3;
        else if(p->gcb4[n][4]==5)
        bc=5;
        else if(p->gcb4[n][4]==6)
        bc=6;
        else if(p->gcb4[n][4]==7)
        bc=7;
        else if(p->gcb4[n][4]==8)
        bc=7;
        else if(p->gcb4[n][4]==21)
        bc=21;
        else if(p->gcb4[n][4]==22)
        bc=21;

        for(q=0;q<p->margin;++q)
        {
            if(p->gcb4[n][3]==1)
            p->flag5[(i-p->imin-1-q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin] = bc;
            else if(p->gcb4[n][3]==2)
            p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1+q)*p->kmax + k-p->kmin] = bc;
            else if(p->gcb4[n][3]==3)
            p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1-q)*p->kmax + k-p->kmin] = bc;
            else if(p->gcb4[n][3]==4)
            p->flag5[(i-p->imin+1+q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin] = bc;
            else if(p->gcb4[n][3]==5)
            p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1-q] = bc;
            else if(p->gcb4[n][3]==6)
            p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1+q] = bc;
        }
    }

    for(n=0;n<p->gcpara1_count;++n)
    {
        i=p->gcpara1[n][0];
        j=p->gcpara1[n][1];
        k=p->gcpara1[n][2];

        p->flag5[IJK]=-10;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin-1-q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=-1;
    }

    for(n=0;n<p->gcpara2_count;++n)
    {
        i=p->gcpara2[n][0];
        j=p->gcpara2[n][1];
        k=p->gcpara2[n][2];

        p->flag5[IJK]=TOPO_FLAG;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1+q)*p->kmax + k-p->kmin]=-2;
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
        i=p->gcpara3[n][0];
        j=p->gcpara3[n][1];
        k=p->gcpara3[n][2];

        p->flag5[IJK]=-10;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1-q)*p->kmax + k-p->kmin]=-3;
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
        i=p->gcpara4[n][0];
        j=p->gcpara4[n][1];
        k=p->gcpara4[n][2];

        p->flag5[IJK]=-10;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin+1+q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=-4;
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
        i=p->gcpara5[n][0];
        j=p->gcpara5[n][1];
        k=p->gcpara5[n][2];

        p->flag5[IJK]=-10;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1-q]=-5;
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
        i=p->gcpara6[n][0];
        j=p->gcpara6[n][1];
        k=p->gcpara6[n][2];

        p->flag5[IJK]=-10;

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1+q]=-6;
    }

    // paraco
    for(n=0;n<p->gcparaco1_count;++n)
    {
        i=p->gcparaco1[n][0];
        j=p->gcparaco1[n][1];
        k=p->gcparaco1[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin-1-q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=-1;
    }

    for(n=0;n<p->gcparaco2_count;++n)
    {
        i=p->gcparaco2[n][0];
        j=p->gcparaco2[n][1];
        k=p->gcparaco2[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1+q)*p->kmax + k-p->kmin]=-2;
    }

    for(n=0;n<p->gcparaco3_count;++n)
    {
        i=p->gcparaco3[n][0];
        j=p->gcparaco3[n][1];
        k=p->gcparaco3[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1-q)*p->kmax + k-p->kmin]=-3;
    }

    for(n=0;n<p->gcparaco4_count;++n)
    {
        i=p->gcparaco4[n][0];
        j=p->gcparaco4[n][1];
        k=p->gcparaco4[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin+1+q)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]=-4;
    }

    for(n=0;n<p->gcparaco5_count;++n)
    {
        i=p->gcparaco5[n][0];
        j=p->gcparaco5[n][1];
        k=p->gcparaco5[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1-q]=-5;
    }

    for(n=0;n<p->gcparaco6_count;++n)
    {
        i=p->gcparaco6[n][0];
        j=p->gcparaco6[n][1];
        k=p->gcparaco6[n][2];

        for(q=0;q<paramargin;++q)
        p->flag5[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1+q]=-6;
    }
}
