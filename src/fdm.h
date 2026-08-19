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

#ifndef FDM_H_
#define FDM_H_

#include"field_header.h"
#include"fieldint1.h"
#include"fieldint2.h"
#include"fieldint3.h"
#include"fieldint4.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"increment.h"
#include"vec.h"
#include"matrix_diag.h"
#include"looping.h"
#include<iostream>
#include<vector>
#if USE_AMREX
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#endif

class lexer;

using namespace std;

class fdm : public increment
{
public:

    fdm(lexer*);

#if USE_AMREX
    // Single MultiFab holding all 47 field components.
    // Declared before the field members so it is constructed first
    // (C++ initialises members in declaration order).
    // Component layout:
    //   0:u  1:v  2:w  3:F  4:G  5:H  6:Fext  7:Gext  8:Hext
    //   9:fbh1  10:fbh2  11:fbh3  12:fbh4  13:fbh5
    //  14:press  15:Fi  16:eddyv  17:L
    //  18:ro  19:visc  20:phi  21:conc  22:test
    //  23:topo  24:solid  25:fb  26:porosity  27:porpart  28:walld
    //  29:nX  30:nY  31:nZ  32:Alpha  33:phasemarker
    //  34:vof  35:vof_nt  36:vof_nb  37:vof_st  38:vof_sb
    //  39:vof_nte  40:vof_ntw  41:vof_nbe  42:vof_nbw
    //  43:vof_ste  44:vof_stw  45:vof_sbe  46:vof_sbw
    amrex::Vector<amrex::MultiFab> m_mf; ///< all fields: 47 components
    amrex::Vector<amrex::MultiFab> m_mf_diag;
#endif

	double gi,gj,gk;

	field1 u,F,Fext;
	field2 v,G,Gext;
	field3 w,H,Hext;
	field4 press;
    field4 Fi;
	field4 eddyv;
	field4 L;
	field4 ro,visc;
	field4 phi;
	field4 conc;
    field4 test;
	field4 topo,solid;
	field4 fb;
	field4 porosity,porpart;
	field4 walld;
	 
	fieldint4 nodeval;
    sliceint4 nodeval2D;
   
    // 6DOF
    field1 fbh1;
    field2 fbh2;
    field3 fbh3;
    field4 fbh4;
    field4 fbh5;
    
    //PLIC
    field4 nX,nY,nZ,Alpha;
    field4 phasemarker;
    field4 vof;
    field4 vof_nt,vof_nb,vof_st,vof_sb; //2D
    field4 vof_nte,vof_ntw,vof_nbe,vof_nbw,vof_ste,vof_stw,vof_sbe,vof_sbw; //3D
    
    // PTF
    slice4 eta,eta_n,depth,WL;
    slice4 Fifsf;
    slice4 K;
    sliceint4 etaloc;
    
    slice1 P;
    slice2 Q;
    
    slice4 bed;
    
	vec rhsvec;

	matrix_diag M;

    double maxF,maxG,maxH;
    double wd_criterion;
	
    #if USE_AMREX
	double t1,t2,t3,t4,t5;
    field1 u0,du0;
    field2 v0,dv0;
    field3 w0,dw0;
    field4 pcorr0,div0,phi0,ro0;
    field1 u1,du1;
    field2 v1,dv1;
    field3 w1,dw1;
    field4 pcorr1,div1,phi1,ro1;
    field1 u2,du2;
    field2 v2,dv2;
    field3 w2,dw2;
    field4 pcorr2,div2,phi2,ro2;
    field4 press0;
    #endif
    field4 grav_pot;

    // Staggered face density rho_face, cell (i,j,k) holding its HIGH face -- the same
    // convention as u/v/w and momentum_FCC3's rox/roy/roz. Filled once per projection by
    // density::update_faces, which also makes the value single-valued across every
    // coarse-fine interface (AMReX MLABecLaplacian does the equivalent internally via
    // averageDownCoeffsToCoarseAmrLevel). Consumers read these instead of recomputing
    // pd->roface(...) from phi per call site, which is what let the coarse and fine sides
    // of a C-F face disagree whenever the smoothed Heaviside band crossed it.
    field1 rofx;
    field2 rofy;
    field3 rofz;

    fieldint4 Mrow;
};

#endif
