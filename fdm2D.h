/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"slice4a.h"
#include"slice5.h"
#include"sliceint1.h"
#include"sliceint2.h"
#include"sliceint4.h"
#include"sliceint5.h"
#include"increment.h"
#include"vec2D.h"
#include"matrix2D.h"
#include"cpt2D.h"
#include"looping.h"
#include"iterators.h"
#include<iostream>

class lexer;

#ifndef FDM2D_H_
#define FDM2D_H_

using namespace std;

class fdm2D : public increment
{
public:

    fdm2D(lexer*);

	double gi,gj,gk;
    
    slice4 eta,eta_n;
    slice1 P,Pn,F;
    slice2 Q,Qn,G;
    slice4 L;
    slice4 ws,wb;
    slice4 press;
    slice4 eddyv,kin,eps;
	slice1 hx;
	slice2 hy;
	slice4 hp,dpx,dpy;
    slice4 test;
    
    slice4 geo,bed,depth;
    slice5 bednode;
    sliceint5 nodeval;
    sliceint4 breaking; 
    slice4 breaking_print;
    
    sliceint4 wet4;
	
	vec2D xvec,rhsvec;

	matrix2D M;
	cpt2D C1,C2,C4;

    double maxF,maxG,maxH,maxK,maxE;
	double inverse,deltax,sigT,Ui,Ua,Uo;
	const double cmu;
	
	double t1,t2,t3,t4,t5;

	int **pvccnode,**tpccnode,*ccedge;
    
    
    void gridsize(lexer*);
};

#endif
