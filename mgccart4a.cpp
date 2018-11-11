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

#include"mgc4a.h"
#include"cart4a.h"
#include"lexer.h"

mgc4a::mgc4a(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
}

mgc4a::~mgc4a()
{
}

void mgc4a::makemgc(lexer* p)
{
	p->Iarray(p->mgc4a,kmax*jmax*imax);

	//make gcdir
	p->gcdirsize4a=1;	
	p->Iarray(p->gcorig4a, p->gcdirsize4a, 6,4);
}

void mgc4a::mgcsetup(lexer* p)
{
    for(i=0;i<imax*jmax*kmax;++i)
	p->mgc4a[i]=0;

	ALOOP
	p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
}


void mgc4a::fillmgc(lexer* p)
{

	int q,n;
	p->gcextra4a=10;
//--------------------------
//WALL1
	QGC4ALOOP
	{
	    i=p->gcb4a[q][0];
		j=p->gcb4a[q][1];
		k=p->gcb4a[q][2];

		if(fabs(p->gcb4a[q][3])==1)
		for(n=0;n<p->margin;++n)
        p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[q][3])==4)
		for(n=0;n<p->margin;++n)
		p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[q][3])==3)
		for(n=0;n<p->margin;++n)
        p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[q][3])==2)
		for(n=0;n<p->margin;++n)
		p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(fabs(p->gcb4a[q][3])==5)
		for(n=0;n<p->margin;++n)
		p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(fabs(p->gcb4a[q][3])==6)
		for(n=0;n<p->margin;++n)
		p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}


//--------------------------
//WALL2
	QGC4ALOOP
	{
        i=p->gcb4a[q][0];
		j=p->gcb4a[q][1];
		k=p->gcb4a[q][2];

		if(fabs(p->gcb4a[q][3])==1)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4a;
			p->gcextra4a++;
        }

		if(fabs(p->gcb4a[q][3])==4)
		for(n=0;n<p->margin;++n)
        if( p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4a;
			p->gcextra4a++;
        }

		if(fabs(p->gcb4a[q][3])==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1
		&& p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]<10)
        {
			p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]=p->gcextra4a;
			p->gcextra4a++;
        }

		if(fabs(p->gcb4a[q][3])==2)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1
		&& p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]<10)
        {
			p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]=p->gcextra4a;
			p->gcextra4a++;
        }

		if(fabs(p->gcb4a[q][3])==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1
		&& p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]<10)
        {
			p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]=p->gcextra4a;
			p->gcextra4a++;
        }

		if(fabs(p->gcb4a[q][3])==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1
		&& p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]<10)
        {
			p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]=p->gcextra4a;
			p->gcextra4a++;
        }
	}
}

void mgc4a::gcdirfill(lexer* p)
{
// GCORIG
    int q,n;
	
	p->Iresize(p->gcorig4a,p->gcdirsize4a, p->gcextra4a, 6, 6, 4, 4); 
	p->gcdirsize4a=p->gcextra4a;
	
	
	for(n=0;n<p->gcdirsize4a;++n)
	for(q=0;q<6;++q)	
	for(qn=0;qn<4;++qn)	
	p->gcorig4a[n][q][qn]=0;	
	
	QGC4LOOP
	{
        i=p->gcb4a[q][0];
		j=p->gcb4a[q][1];
		k=p->gcb4a[q][2];

		if(p->gcb4a[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][0][di]=1;
        }

		if(p->gcb4a[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1)
        {
			di = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]-10][3][di]=1;
        }

		if(p->gcb4a[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]-10][2][dj]=1;
	    }

		if(p->gcb4a[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1)
        {
			dj = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]-10][1][dj]=1;
        }

		if(p->gcb4a[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1)
        {
			dk = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]-10][4][dk]=1;
        }

		if(p->gcb4a[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1)
        {
			dk = (n+1);
			p->gcorig4a[p->mgc4a[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]-10][5][dk]=1;
        }
	}


}



void mgc4a::fillgcb(lexer *p)
{
    int q;
    
    QGC4ALOOP
	{
	for(n=0;n<5;++n)
	p->gcb4a[q][n]=p->gcb4[q][n];
    
    p->gcd4a[q]=p->gcd4[q];
	}
}


