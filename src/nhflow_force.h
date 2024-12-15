/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_FORCE_H_
#define NHFLOW_FORCE_H_

#include"fieldint5.h"
#include"field5.h"
#include"increment.h"
#include"vtp3D.h"
#include<iostream>
#include<fstream>
#include<vector>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class nhflow_force : virtual public increment, private vtp3D
{

public:
	nhflow_force(lexer*,fdm_nhf*,ghostcell*,int);
	virtual ~nhflow_force();
	void start(lexer*,fdm_nhf*,ghostcell*);
    void ini(lexer*,fdm_nhf*,ghostcell*);

private:
	void triangulation(lexer*, fdm_nhf*, ghostcell*);
	void reconstruct(lexer*, fdm_nhf*);
	void addpoint(lexer*,fdm_nhf*,int,int);
    
    void allocate(lexer*,fdm_nhf*,ghostcell*);
    void deallocate(lexer*,fdm_nhf*,ghostcell*);

    std::vector<char> buffer;
    int m = 0;

    int *vertice,*nodeflag;
    double *eta;
	
	int **tri, **facet, *confac, *numfac,*numpt;
	double **ccpt, **pt, *ls;
	double   dV1,dV2,C1,C2,mi;
	int numtri,numvert, numtri_mem, numvert_mem;
	int count,countM,n,nn,q;
	int ccptcount,facount,check;
	int polygon_sum,polygon_num,point_num;
	const double zero,interfac;
    double epsi;
	
	
    void force_calc(lexer*,fdm_nhf*,ghostcell*);
    
    
	void print_force(lexer*,fdm_nhf*,ghostcell*);
    void print_ini(lexer*,fdm_nhf*,ghostcell*);
    void print_vtp(lexer*,fdm_nhf*,ghostcell*);
    void pvtp(lexer*,fdm_nhf*,ghostcell*);
    void name_iter(lexer*,fdm_nhf*,ghostcell*);
    void piecename(lexer*,fdm_nhf*,ghostcell*,int);

    char name[100],pname[100];
    int iin,offset[100];
    float ffn;
	int forceprintcount;
    int gcval_press;
    
    // force
    double Fx,Fy,Fz;
    double A_tot,A;
    
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double xc,yc,zc;
    double nx,ny,nz,norm;
    double nxs,nys,nzs;
    double uval,vval,wval,pval,viscosity,density,phival;
    double etaval,hspval;
    double du,dv,dw;
    double at,bt,ct,st;
    
    ofstream fout;
    
    double xs,xe,ys,ye,zs,ze;
    double xm,ym,zm;
	int is,ie,js,je,ks,ke;
    const int ID;
	
    double ux,vy,wz,vel,pressure;
    double xloc,yloc,zloc;
	double xlocvel,ylocvel,zlocvel;
    double sgnx,sgny,sgnz;
    double Ax;
    double Ay;
    double Az;

    double xp1,xp2,yp1,yp2,zp1,zp2;
};

#endif


