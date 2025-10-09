/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

void ghostcell::gcslflagx(lexer* p, int *flag)
{
    count=0;
    for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];

        isend1[count]=flag[(i-p->imin)*p->jmax + (j-p->jmin)];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];

        isend2[count]=flag[(i-p->imin)*p->jmax + (j-p->jmin)];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];

        isend3[count]=flag[(i-p->imin)*p->jmax + (j-p->jmin)];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];

        isend4[count]=flag[(i-p->imin)*p->jmax + (j-p->jmin)];
        ++count;
    }

    Sendrecv_int(p->gcslpara1_count, p->gcslpara2_count, p->gcslpara3_count, p->gcslpara4_count, 0, 0);

//  Unpack

    count=0;
    for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];
    k=p->gcslpara1[n][2];

        flag[(i-p->imin-1)*p->jmax + (j-p->jmin)]=irecv1[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
    k=p->gcslpara2[n][2];

        flag[(i-p->imin)*p->jmax + (j-p->jmin+1)]=irecv2[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];
    k=p->gcslpara3[n][2];

        flag[(i-p->imin)*p->jmax + (j-p->jmin-1)]=irecv3[count];
        ++count;
    }

    count=0;
    for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];
    k=p->gcslpara4[n][2];

        flag[(i-p->imin+1)*p->jmax + (j-p->jmin)]=irecv4[count];
        ++count;
    }
}
