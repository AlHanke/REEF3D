/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef FORCE_H_
#define FORCE_H_

#include"fieldint5.h"
#include"field5.h"
#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;

using namespace std;

class force : public increment
{

public:
    force(lexer*,fdm*,ghostcell*,int);
    virtual ~force() = default;
    virtual void start(lexer*,fdm*,ghostcell*);
    virtual void ini(lexer*,fdm*,ghostcell*);

private:
    void triangulation(lexer*,fdm*);
    void reconstruct(lexer*,fdm*);
    void addpoint(lexer*,int,int,int);
    void force_calc(lexer*,fdm*,ghostcell*);
    
    void print_force(lexer*);
    void print_ini(lexer*);
    void print_vtp(lexer*,fdm*);
    void pvtp(int);

    int numtri,numvert;
    int n,q;
    int ccptcount,facount;
    int is,ie,js,je,ks,ke;
    int gcval_press;

    const int ID;
    const double interfac;
    const double threshold;

    double xs,xe,ys,ye,zs,ze;
    double xm,ym,zm;

    // force
    double Fx,Fy,Fz;
    double A_tot;
    
    int **tri, **facet, *confac, *numfac,*numpt;
    double **ccpt, **pt, *ls;
    
    fieldint5 vertice, nodeflag;
    field5 eta;

    char name[100],pname[100];
    int iin,offset[100];
    float ffn;
    
    ofstream fout;
};

#endif
