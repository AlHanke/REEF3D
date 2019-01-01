/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

int ghostcell::column_pt1_count(lexer* p, fdm* a)
{
	n=0;
    ULOOP
	++n;
	
	
	GGC1LOOP
    {
    i=p->gcb1[g][0];
    j=p->gcb1[g][1];
    k=p->gcb1[g][2];

    for(q=0;q<margin;++q)
    ++n;
    
    }

	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)		
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
    
    //cout<<p->mpirank<<"Cpt1_n: "<<n<<"  p->veclength: "<<p->veclength<<endl;
    
    return n;
}

int ghostcell::column_pt2_count(lexer* p, fdm* a)
{
	n=0;
    VLOOP
	++n;

	
	//int cellnum1=n;
	GGC2LOOP
    {
    i=p->gcb2[g][0];
    j=p->gcb2[g][1];
    k=p->gcb2[g][2];

    for(q=0;q<margin;++q)
    ++n;
    }
	
	//cout<<p->mpirank<<" gcbcount2: "<<p->gcb4_count*3<<" N2: "<<n-cellnum1<<endl;
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
	
    return n;
}

int ghostcell::column_pt3_count(lexer* p, fdm* a)
{
	n=0;
    WLOOP
	++n;
	
	GGC3LOOP
    {
    i=p->gcb3[g][0];
    j=p->gcb3[g][1];
    k=p->gcb3[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;

	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
    
    return n;
}

int ghostcell::column_pt4_count(lexer* p, fdm* a)
{
	n=0;
    LOOP
	++n;

	
	//int cellnum1=n;
	GGC4LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	//cout<<p->mpirank<<" gcbcount4: "<<p->gcb4_count*3<<" N4: "<<n<<" cellnum: "<<p->cellnum<<endl;
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
    
    return n;
}

int ghostcell::column_pt4a_count(lexer* p, fdm* a)
{
	n=0;
    ALOOP
	++n;
	
	//int cellnum1=n;
	GGC4ALOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		++n;
	}
    
    return n;
}

int ghostcell::column_pt6_count(lexer* p, fdm* a)
{
	n=0;
    BASELOOP
	++n;

	
	//int cellnum1=n;
	GGC6LOOP
    {
    i=p->gcb4[g][0];
    j=p->gcb4[g][1];
    k=p->gcb4[g][2];

        for(q=0;q<margin;++q)
        ++n;
    }
	
	//cout<<p->mpirank<<" gcbcount4: "<<p->gcb4_count*3<<" N4: "<<n-cellnum1<<endl;
	
	for(g=0;g<p->gcpara1_count;++g)
    {
    i=p->gcpara1[g][0];
    j=p->gcpara1[g][1];
    k=p->gcpara1[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara1[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara2_count;++g)
    {
    i=p->gcpara2[g][0];
    j=p->gcpara2[g][1];
    k=p->gcpara2[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara2[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara3_count;++g)
    {
    i=p->gcpara3[g][0];
    j=p->gcpara3[g][1];
    k=p->gcpara3[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara3[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara4_count;++g)
    {
    i=p->gcpara4[g][0];
    j=p->gcpara4[g][1];
    k=p->gcpara4[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara4[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara5_count;++g)
    {
    i=p->gcpara5[g][0];
    j=p->gcpara5[g][1];
    k=p->gcpara5[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara5[g][3]==1)
		++n;
	}
	
	for(g=0;g<p->gcpara6_count;++g)
    {
    i=p->gcpara6[g][0];
    j=p->gcpara6[g][1];
    k=p->gcpara6[g][2];
	
		for(q=0;q<margin;++q)
		if(p->gcpara6[g][3]==1)
		++n;
	}
    
    return n;
}