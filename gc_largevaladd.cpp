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
#include"field.h"
#include<math.h>

void ghostcell::largevaladd(field &f,double dist,int gcv, int bc, int cs)
{
	int qstart=1;
	
	margin=3;

	if(cs==1)
	for(q=qstart;q<margin;++q)
	f(i-q-1,j,k)=pow(-1.0,q+1)*double(q+1)*1.0e20;

	if(cs==2)
	for(q=qstart;q<margin;++q)
	f(i,j+q+1,k)=pow(-1.0,q+1)*double(q+1)*1.0e20;

	if(cs==3)
	for(q=qstart;q<margin;++q)
	f(i,j-q-1,k)=-pow(-1.0,q+1)*double(q+1)*1.0e20;

	if(cs==4)
	for(q=qstart;q<margin;++q)
	f(i+q+1,j,k)=-pow(-1.0,q+1)*double(q+1)*1.0e20;

	if(cs==5)
	for(q=qstart;q<margin;++q)
	f(i,j,k-q-1)=-pow(-1.0,q+1)*double(q+1)*1.0e20;

	if(cs==6)
	for(q=qstart;q<margin;++q)
	f(i,j,k+q+1)=pow(-1.0,q+1)*double(q+1)*1.0e20;
}
