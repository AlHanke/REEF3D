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

#include"mgcslice2.h"
#include"lexer.h"
#include"fdm.h"

void mgcslice2::make_ggc(lexer* p)
{
	p->ggcslsize2=1;
	p->Iarray(p->ggcsl2,p->ggcslsize2,2);
}

void mgcslice2::fill_ggc(lexer* p)
{
	int q,qq,n,nn,a;
	int check;
	
    int *ggcslmem2;
	p->Iarray(ggcslmem2,imax*jmax);

//--------------------------
//WALL1

	GCSL2LOOP
	{
        i=p->gcbsl2[n][0];
		j=p->gcbsl2[n][1];

		if(p->gcbsl2[n][3]==1)
		for(q=0;q<p->margin;++q)
            ggcslmem2[(i-imin-q-1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl2[n][3]==4)
		for(q=0;q<p->margin;++q)
		    ggcslmem2[(i-imin+q+1)*jmax + (j-jmin)]+=1;

		if(p->gcbsl2[n][3]==3)
		for(q=0;q<p->margin;++q)
            ggcslmem2[(i-imin)*jmax + (j-jmin-q-1)]+=1;

		if(p->gcbsl2[n][3]==2)
		for(q=0;q<p->margin;++q)
		    ggcslmem2[(i-imin)*jmax + (j-jmin+q+1)]+=1;

	}

// count entries
	p->ggcslcount2=0;
	a=0;
	for(i=0;i<imax;++i)
	for(j=0;j<jmax;++j)
	{
        if(ggcslmem2[a]>1)
        ++p->ggcslcount2;

	++a;
	}
	
	p->Iresize(p->ggcsl2,p->ggcslsize2,p->ggcslcount2*p->margin, 2, 2);
	p->ggcslsize2=p->ggcslcount2*p->margin;

//--------------------------
//WALL2
	n=0;
	QQGCSL2LOOP
	{
        i=p->gcbsl2[qq][0];
		j=p->gcbsl2[qq][1];

		if(p->gcbsl2[qq][3]==1)
		for(q=0;q<p->margin;++q)
		if(ggcslmem2[(i-imin-q-1)*jmax + (j-jmin)]>1)
        {
            if(ggcslmem2[(i-imin-q-1)*jmax + (j-jmin)]<10)
            {
                ggcslmem2[(i-imin-q-1)*jmax + (j-jmin)]=n+10;
			 p->ggcsl2[n][0]=i-q-1;
			 p->ggcsl2[n][1]=j;
			 ++n;
            }
        }

        if(p->gcbsl2[qq][3]==4)
        for(q=0;q<p->margin;++q)
		if(ggcslmem2[(i-imin+q+1)*jmax + (j-jmin)]>1)
        {
            if(ggcslmem2[(i-imin+q+1)*jmax + (j-jmin)]<10)
            {
                ggcslmem2[(i-imin+q+1)*jmax + (j-jmin)]=n+10;
			p->ggcsl2[n][0]=i+q+1;
			p->ggcsl2[n][1]=j;
			++n;
            }
        }

        if(p->gcbsl2[qq][3]==3)
        for(q=0;q<p->margin;++q)
		if(ggcslmem2[(i-imin)*jmax + (j-jmin-q-1)]>1)
        {
            if(ggcslmem2[(i-imin)*jmax + (j-jmin-q-1)]<10)
            {
                ggcslmem2[(i-imin)*jmax + (j-jmin-q-1)]=n+10;
			p->ggcsl2[n][0]=i;
			p->ggcsl2[n][1]=j-q-1;
			++n;
            }
        }

		if(p->gcbsl2[qq][3]==2)
		for(q=0;q<p->margin;++q)
		if(ggcslmem2[(i-imin)*jmax + (j-jmin+q+1)]>1)
        {
            if(ggcslmem2[(i-imin)*jmax + (j-jmin+q+1)]<10)
            {
                ggcslmem2[(i-imin)*jmax + (j-jmin+q+1)]=n+10;
			p->ggcsl2[n][0]=i;
			p->ggcsl2[n][1]=j+q+1;
			++n;
            }
        }
	}
	p->ggcslcount2=n;


	p->del_Iarray(ggcslmem2,imax*jmax);
}


