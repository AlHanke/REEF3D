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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"

void ioflow_f::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2;

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][5]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
    
    if(p->I10==1)
    velini(p,a,pgc);
    
    // 4a ---------------
	
	count1=0;
    count2=0;
    GC4ALOOP
    {
        if(p->gcb4a[n][4]==1)
        {
        p->gcin4a[count1][0]=p->gcb4a[n][0];
        p->gcin4a[count1][1]=p->gcb4a[n][1];
        p->gcin4a[count1][2]=p->gcb4a[n][2];
        p->gcin4a[count1][3]=p->gcb4a[n][3];
        p->gcin4a[count1][5]=p->gcb4a[n][5];
        ++count1;
        }

        if(p->gcb4a[n][4]==2)
        {
        p->gcout4a[count2][0]=p->gcb4a[n][0];
        p->gcout4a[count2][1]=p->gcb4a[n][1];
        p->gcout4a[count2][2]=p->gcb4a[n][2];
        p->gcout4a[count2][3]=p->gcb4a[n][3];
        p->gcout4a[count2][5]=p->gcb4a[n][5];
        ++count2;
        }
    }

    p->gcin4a_count=count1;
    p->gcout4a_count=count2;
    

}

void ioflow_f::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, discrete *pdisc, reini *preini, ioflow *pflow)
{

	p->del_Darray(walldin, walldin_size);
	p->del_Darray(walldout, walldout_size);
	
	walldin_size=p->gcin_count;
	walldout_size=p->gcout_count;
	
	p->Darray(walldin, walldin_size);
    p->Darray(walldout, walldout_size);


    for(n=0;n<p->gcin_count;++n)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    walldin[n] = a->walld(i,j,k);
    }

    for(n=0;n<p->gcout_count;++n)
    {
    i=p->gcout[n][0];
    j=p->gcout[n][1];
    k=p->gcout[n][2];

    walldout[n] = a->walld(i,j,k);
    }
}


void ioflow_f::periodic(field& b, lexer *p)
{
	if(p->B210==1 && p->count>0)
	{
        for(int n=0;n<p->gcb4_count;++n)
        if(p->gcb4[n][4]==1)
        {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        b(i-1,j,k)=b(p->knox-1,j,k);
        b(i-2,j,k)=b(p->knox-1,j,k);
        b(i-3,j,k)=b(p->knox-1,j,k);
        }
	}
}

void ioflow_f::iogcb_update(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2;

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][5]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
}

void ioflow_f::veltimesave(lexer *p, fdm *a, ghostcell *pgc)
{
    pvrans->veltimesave(p,a,pgc);
}