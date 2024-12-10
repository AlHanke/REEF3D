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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef PART_H_
#define PART_H_

enum {EMPTY=-1,PASSIVE=1,ACTIVE=10,MOVING=20};


#include"increment.h"

class lexer;
class ghostcell;
class slice;

class part : public increment
{
public:
    part(lexer*, ghostcell *);
    ~part();
    
// functions
    // add
    void add(lexer*,ghostcell*,double,double,double,double,double);
    
    void resize(lexer*,int);

    // remove
    void remove(int);
    
    // parallel
    void xchange(lexer*, ghostcell*,slice&,int);
    
// data arrays
    double *U,*V,*W;
    double *URK1,*VRK1,*WRK1;
    
    double *X,*Y,*Z;
    double *XRK1,*YRK1,*ZRK1;
    
    double *Uf,*Vf,*Wf;

    double *D,*RO;
    
    double* Test;
    double* Test2;
    
    int *Empty,*Flag;
    
    double d50,rhosed;
    double ParcelFactor;
    
// iterators
    int index; // replace loopindex
    int index_empty; //
    
    int n,q;
    
private:
    void xchange_count(lexer*, ghostcell*,int);
    void xchange_sendid(lexer*, ghostcell*,int);
    void xchange_fill(lexer*, ghostcell*, int, double*);
    void xchange_fill_flag(lexer*, ghostcell*,int);
    void xchange_fillback(lexer*, ghostcell*, double*);
    void xchange_fillback_flag(lexer*, ghostcell*,slice&,int);
    void xchange_resize(lexer*, ghostcell*);

    int sendnum[6];
    int recvnum[6];
    int sendcount[6];
    int recvcount[6];
    int** sendid;
    double** send;
    double** recv;
    int capacity_para;
    int maxnum;
    int index_empty0;
    int capacity; // length of allocated array
};

#endif
