/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_convection.h"
#include"slice1.h"
#include"slice2.h"
#include"increment.h"

class nhflow_flux;
class nhflow_reconstruct;
class nhflow_signal_speed;
class patchBC_interface;
class ghostcell;

#ifndef NHFLOW_HLL_H_
#define NHFLOW_HLL_H_

using namespace std;

class nhflow_HLL : public nhflow_convection, public increment
{

public:

	nhflow_HLL (lexer*,ghostcell*,patchBC_interface*);
	virtual ~nhflow_HLL();

    virtual void start(lexer*, fdm_nhf*, double*, int, double*, double*,double*, slice&);
    virtual void precalc(lexer*, fdm_nhf*, double*, int, double*, double*,double*, slice&);

private:
    double aij(lexer*, fdm_nhf*, double*, int, double*, double*, double*, double*, double*, double*);
    
    double aij_U(lexer*, fdm_nhf*, double*, int, double*, double*, double*, double*, double*, double*);
    double aij_V(lexer*, fdm_nhf*, double*, int, double*, double*, double*, double*, double*, double*);
    double aij_W(lexer*, fdm_nhf*, double*, int, double*, double*, double*, double*, double*, double*);
    

    double *Fs,*Fn,*Fe,*Fw,*Fz;
    double *Us,*Un,*Ue,*Uw,*Ub,*Ut;
    double *Vs,*Vn,*Ve,*Vw,*Vb,*Vt;

    slice1 ETAs,ETAn;
    slice2 ETAe,ETAw;
    slice1 Ds,Dn;
    slice2 De,Dw;
    slice1 Ss,Sn;
    slice2 Se,Sw;
    
    

	double dx,dy,dz;
	double udir,vdir,wdir;
	double L;

    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    ghostcell *pgc;
    nhflow_flux *pflux;
    patchBC_interface *pBC;
    nhflow_reconstruct *precon;
    nhflow_signal_speed *pss;
};

#endif
