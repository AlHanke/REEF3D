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

#ifndef LEXER_H_
#define LEXER_H_

#include "ArrayWrapper2D.h"
#include "ArrayWrapper3D.h"
#include "control.h"
#include "coordinates.h"
#if USE_AMREX
    #include "grid_amrex.h"
#else
    #include "grid.h"
#endif
#include "increment.h"
#include "interpolation.h"
#include "looping.h"
#include "position.h"
#include "resize.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>

class ghostcell;
class weno_nug_func;

using namespace std;

class lexer
    : public virtual resize_class,
      public control,
      public coordinates,
      #if USE_AMREX
          public grid_amrex,
      #else
          public grid,
      #endif
      public interpolation,
      public position
{
public:
    lexer();
    ~lexer() override = default;

//-----functions------------------
//---- setup
    void lexer_read(ghostcell*);
    void flagini();
    void gridini(ghostcell*);
    void makeflag(int*);

    void read_grid();
    void control_calc();
    void gridsize();
    void vecsize(ghostcell*);
    void vellast();
    void lexer_ini();
    int conv(double);

    // 2D
    void grid2Dsize();
    void flagini2D();
    void gridini2D();

//-----data-----------------------
    //REEF3D

    int pointnum,cellnum,tpcellnum;
    int cellnumtot;
    int N4,N4_row;
    int N7,N7_row;
    int *flag1,*flag2,*flag3,*flag4,*flag5,*flag7;

    // flag
    std::unique_ptr<int[]> flag_solid, flag_topo;
    std::unique_ptr<double[]> data;
    std::unique_ptr<double[]> topobed, solidbed;
    double *bed,*WL;
    int *wet;
    int *deep;
    int gcbextra;
    int solidread,toporead,topoforcing;
    int cms_flag;

    //GHOSTCELL
    int **gcb1,**gcb2,**gcb3,**gcb4;

    int gcdf1_count,gcdf2_count,gcdf3_count,gcdf4_count;
    std::vector<std::array<int, 4>> gcdf1, gcdf2, gcdf3;
    std::vector<std::array<int, 5>> gcdf4;
    int gcsldfeta4_count,gcsldfbed4_count;
    std::vector<std::array<int,3>> gcsldfeta4, gcsldfbed4;

    int gcb1_count,gcb2_count,gcb3_count,gcb4_count;
    int gcpara_sum;
    int solid_gcb_est, topo_gcb_est, solid_gcbextra_est, topo_gcbextra_est, tot_gcbextra_est;
    int bcside1,bcside2,bcside3,bcside4,bcside5,bcside6;

    // serial periodic BC
    int periodic1,periodic2,periodic3;
    int periodicX1,periodicX2,periodicX3,periodicX4,periodicX5,periodicX6;

    int **dgc1,**dgc2,**dgc3,**dgc4;
    int dgc1_count,dgc2_count,dgc3_count,dgc4_count;

    // PARALLEL
    std::vector<std::array<int,7>> gcpara1;
    std::vector<std::array<int,7>> gcpara2;
    std::vector<std::array<int,7>> gcpara3;
    std::vector<std::array<int,7>> gcpara4;
    std::vector<std::array<int,7>> gcpara5;
    std::vector<std::array<int,7>> gcpara6;

    std::vector<std::array<int,3>> gcparaco1;
    std::vector<std::array<int,3>> gcparaco2;
    std::vector<std::array<int,3>> gcparaco3;
    std::vector<std::array<int,3>> gcparaco4;
    std::vector<std::array<int,3>> gcparaco5;
    std::vector<std::array<int,3>> gcparaco6;

    std::array<std::vector<std::array<int,3>>, 4> gcx7;
    std::array<int, 4> gcx7_count;
    std::array<std::vector<std::array<int,3>>, 4> gcxco7;
    std::array<int, 4> gcxco7_count;

    int gcpara1_count, gcpara2_count, gcpara3_count, gcpara4_count, gcpara5_count, gcpara6_count;
    int gcparaco1_count, gcparaco2_count, gcparaco3_count, gcparaco4_count, gcparaco5_count, gcparaco6_count;
    int gcslpara1_count, gcslpara2_count, gcslpara3_count, gcslpara4_count;
    int gcslparaco1_count, gcslparaco2_count, gcslparaco3_count, gcslparaco4_count;
    int nb1,nb2,nb3,nb4,nb5,nb6;
    int mx,my,mz;
    int mpi_size;

    int ulast,vlast,wlast,flast;

    // Solver
    std::vector<int> range_col4, range_row4, range_col7, range_row7;
    int sizeM4[2] = {0,0};

    // SMO
    int veclength;

    //SLICE
    int *flagslice1,*flagslice2,*flagslice4;

    int vec2Dlength;

    int pointnum2D,cellnum2D,cellnumtot2D,polygon_sum;

    // SLICE ghostcell
    int gcbsl1_count,gcbsl2_count,gcbsl4_count;
    int gcslin_count,gcslout_count;
    int gcslawa1_count,gcslawa2_count;
    int **gcbsl1,**gcbsl2,**gcbsl4;
    int **gcslin, **gcslout;
    int **gcslawa1, **gcslawa2;

    std::vector<std::array<int,3>> dgcsl1, dgcsl2, dgcsl4;
    int dgcsl1_count,dgcsl2_count,dgcsl4_count;

    // SLICE parallel
    std::vector<std::array<int,2>> gcslpara1;
    std::vector<std::array<int,2>> gcslpara2;
    std::vector<std::array<int,2>> gcslpara3;
    std::vector<std::array<int,2>> gcslpara4;

    std::vector<std::array<int,2>> gcslparaco1;
    std::vector<std::array<int,2>> gcslparaco2;
    std::vector<std::array<int,2>> gcslparaco3;
    std::vector<std::array<int,2>> gcslparaco4;

    // flow parameters
    static constexpr double cmu = 0.09;
    double Ui,Ua,Uo;
    double Hi,Ho;
    double Qi,Qo;

    // 6DOF
    double ufb,vfb,wfb;
    double pfb,qfb,rfb;
    double ufbi,vfbi,wfbi;
    double pfbi,qfbi,rfbi;
    double xg,yg,zg;
    double phi_fb,theta_fb,psi_fb;
    double ufbmax, vfbmax, wfbmax;
    int mooring_count, net_count;

    // FSI
    int FSI_count;

    // time + iterations
    int inneriter,count,solveriter,preconiter,count_statestart;
    int solver_status,solver_error;
    int sediter;
    double final_res;
    double dt,dt_old,simtime,viscmax;
    double mindt,maxdt;
    double umax,vmax,wmax,epsmax,kinmax,pressmin,pressmax,omegamax;
    double presstime,reinitime,itertime;
    double wavecalctime;
    double meantime,totaltime;
    double gcmeantime,gctotaltime;
    double Xmeantime,Xtotaltime;
    double susptime;
    double gctime, xtime;
    double volume1,volume2,volume3;
    double dtsed,sedtime,slidecells;
    double printtime,sedprinttime,fsfprinttime,fsfsedprinttime,probeprinttime,stateprinttime;
    double wavetime;

    // solver watch
    int uiter,viter,witer;
    int kiniter,epsiter;
    int poissoniter, laplaceiter;
    int lsmiter;
    int suspiter,topoiter;
    int heatiter,concentrationiter;
    int printcount, printcount_sixdof;
    double utime,vtime,wtime;
    double recontime,fsftime;
    double dftime;
    double kintime,epstime;
    double poissontime, laplacetime, matrixtime, ptime;
    double sftime,fbtime;
    double fbdt,fbmax;
    double sfmax;
    double lsmtime,concentrationtime;
    double printouttime;
    double phimean,phiout,phiin;
    double fsfin,fsfout;
    double pressgage;

    // wave coefficients
    double wT,wV,wH,wA,wL,wd,ww,wk,wC;
    double wHs,wAs,wwp,ww_s,ww_e,wTp,wLp;
    int wN;
    double wts,wte;

    // free surface
    double psi;

// PARALELL
    int mpirank;

    weno_nug_func *wenofunc;

// sigma coordinate
    double *sig;
    double *sigx,*sigy,*sigz,*sigt;
    double *sigxx;
};

#include "ArrayWrapper2D_imp.h"
#include "ArrayWrapper3D_imp.h"

#endif
