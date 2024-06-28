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

#ifndef LEXER_H_
#define LEXER_H_

#include<iostream>
#include<cstdlib>
#include<iomanip>
#include<math.h>
#include"resize.h"
#include"increment.h"
#include"position.h"
#include"interpolation.h"
#include<fstream>
#include"looping.h"
#include<vector>

class weno_nug_func;
class ghostcell;

using namespace std;

class lexer : public increment, public resize_class, public position, public interpolation
{
public:

	lexer();
	virtual ~lexer();

//-----functions------------------
//---- setup
    void lexer_read(ghostcell*);
    void flagini();
	void gridini(ghostcell*);
    void gcd_ini(ghostcell*);
    void gridini_patchBC();
    void makeflag(int*);
    void sigma_coord_ini();
	
	void read_grid();
	void read_control();
	void control_calc();
	void ini_default();
	void assign_margin();
	void ctrlsend();
	void ctrlrecv();
	int maxparacount();
	void gridsize();
	void vecsize(ghostcell*);
    void gcbextra_est(ghostcell*);
	void vellast();
	void indices_minmax();
	void lexer_ini();
    void lexer_gridspacing(ghostcell*);
	void parse();
	void fieldlogic();
	int conv(double);
    
    // 2D
    void grid2Dsize();
    void flagini2D();
	void gridini2D();


//-----data-----------------------

	//REEF3D
	double dx,dy,dz;
    double *xpoint,*ypoint,*zpoint;
    double *xnode,*ynode,*znode;
    
    
	int imin,imax,jmin,jmax,kmin,kmax;
    int kmaxF;
	int pointnum,cellnum,tpcellnum;
	int cellnum1,cellnum2,cellnum3;
    int pointnumtot,cellnumtot;
    int N4,N4_row,N4_col; //time scheme
    int N7,N7_row,N7_col;
	double originx,originy,originz;
    double endx,endy,endz;
	double global_xmin,global_ymin,global_zmin;
	double global_xmax,global_ymax,global_zmax;
	int origin_i, origin_j, origin_k;
	int gknox,gknoy,gknoz;
	int surf_tot;
	int *flag1,*flag2,*flag3,*flag4,*flag5,*flag7,*flag;
    int *flagsf1,*flagsf2,*flagsf3,*flagsf4;
    int *BC;
	int*mgflag;
    double *flag_solid,*flag_topo;
    double *data;
	double *topobed,*solidbed,*bed,*depth;
    int *wet,*wet_n;
    int *deep;
	int *tpflag,*ndbaseflag;
	int *mgc1,*mgc2,*mgc3,*mgc4,*mgc4a,*mgc6;
	int ***gcorig1,***gcorig2,***gcorig3,***gcorig4,***gcorig4a,***gcorig6;
	int gcdirsize1,gcdirsize2,gcdirsize3,gcdirsize4,gcdirsize4a,gcdirsize6;
	int i_dir,j_dir,k_dir;
	double x_dir,y_dir,z_dir;
    int gcbextra;
    int solidread,toporead,porousread;


    //GHOSTCELL
	int **gcb1,**gcb2,**gcb3,**gcb4,**gcb4a,*gcb6;
	int **gcin, **gcout, **gcpress,**gcin6, **gcout6;
	int **gcin4a, **gcout4a;
	double *gcd1,*gcd2,*gcd3,*gcd4,*gcd4a;
	double **gcn;
	int *gcside4;
	int gcside4_size;
	int gcextra1,gcextra2,gcextra3,gcextra4,gcextra4a,gcextra6;
    
    int gcdf1_count,gcdf2_count,gcdf3_count,gcdf4_count;
    int **gcdf1,**gcdf2,**gcdf3,**gcdf4;

	int **dgc1,**dgc2,**dgc3,**dgc4;
	int dgc1_count,dgc2_count,dgc3_count,dgc4_count;

	int gcwall_count, gcin_count, gcout_count, gcpress_count, gcfsf_count, gcbed_count;
    int gcin6_count, gcout6_count;
	int gcin4a_count, gcout4a_count;
	int gcb1_count,gcb2_count,gcb3_count,gcb4_count,gcb4a_count;
	int gcpara_sum, gcparaco_sum;
	int gcb_fix,gcb_solid,gcb_topo,gcb_fb, solid_gcb_est, topo_gcb_est, solid_gcbextra_est, topo_gcbextra_est, tot_gcbextra_est;
	int gcb_sediment_est, gcb_floating_est;
    int bcside1,bcside2,bcside3,bcside4,bcside5,bcside6;
    
    // serial periodic BC
    int periodic1,periodic2,periodic3;
    int periodicX1,periodicX2,periodicX3,periodicX4,periodicX5,periodicX6;
    
    int **gc4periodic;
    int **gc4aperiodic;
    int *gc4periodic_count;
    int *gc4aperiodic_count;
    int gc4periodic_maxcount;
    
	// PARALLEL
	int** gcpara1;
	int** gcpara2;
	int** gcpara3;
	int** gcpara4;
	int** gcpara5;
	int** gcpara6;

	int** gcparaco1;
	int** gcparaco2;
	int** gcparaco3;
	int** gcparaco4;
	int** gcparaco5;
	int** gcparaco6;
    
    int*** gcx7;
    int* gcx7_count;
    int*** gcxco7;
    int* gcxco7_count;

	int gcpara1_count, gcpara2_count, gcpara3_count, gcpara4_count, gcpara5_count, gcpara6_count;
	int gcparaco1_count, gcparaco2_count, gcparaco3_count, gcparaco4_count, gcparaco5_count, gcparaco6_count;
    int gcslpara1_count, gcslpara2_count, gcslpara3_count, gcslpara4_count;
    int gcslparaco1_count, gcslparaco2_count, gcslparaco3_count, gcslparaco4_count;
	int maxpara;
	int nb1,nb2,nb3,nb4,nb5,nb6;
	int mx,my,mz;
    int mi,mj,mk;
	int mpi_edgenum,mpi_nodes,mpi_size;
	int *mpi_index, *mpi_edges;
	
	int ulast,vlast,wlast,flast,ulastsflow;
	int velcorr;
	int* ictrl;
	double* dctrl;
	int ctrlsize;
	int stencil;	

	// Solver
	int *colnum;
    int *range_col4,*range_row4,*range_col7,*range_row7;
	int *sizeM1,*sizeM2,*sizeM3,*sizeM4,*sizeM4a,*sizeM6,*sizeM9;
    int *sizeS1,*sizeS2,*sizeS4; 
	int mglevel_max,*MGL;

	// SMO
	int veclength;
    int C4_size,C4a_size,C6_size;
    int C1_2D_size,C2_2D_size,C4_2D_size;
    int M_size,M_2D_size;
    
    //SLICE
    int *flagslice1,*flagslice2,*flagslice4,*tpflagslice;
    int *flagfsf;
    int *mgcsl1,*mgcsl2,*mgcsl3,*mgcsl4,*mgcsl4a;
    int ***gcslorig1,***gcslorig2,***gcslorig3,***gcslorig4,***gcslorig4a;
	int gcsldirsize1,gcsldirsize2,gcsldirsize3,gcsldirsize4,gcsldirsize4a;
    
    int slicenum,vec2Dlength;
    
    int pointnum2D,cellnum2D,cellnumtot2D,polygon_sum;
    
    // SLICE ghostcell
    int gcbsl1_count,gcbsl2_count,gcbsl3_count,gcbsl4_count,gcbsl4a_count;
    int gcslin_count,gcslout_count;
    int gcslawa1_count,gcslawa2_count;
    int **gcbsl1,**gcbsl2,**gcbsl3,**gcbsl4,**gcbsl4a;
	int **gcslin, **gcslout;
    int **gcslawa1, **gcslawa2;
	double *gcdsl1,*gcdsl2,*gcdsl3,*gcdsl4,*gcdsl4a;


    int gcsl_extra1,gcsl_extra2,gcsl_extra3,gcsl_extra4,gcsl_extra4a;

	int **dgcsl1,**dgcsl2,**dgcsl3,**dgcsl4;
	int dgcsl1_count,dgcsl2_count,dgcsl3_count,dgcsl4_count;
    
    int **ggcsl1,**ggcsl2,**ggcsl3,**ggcsl4,**ggcsl4a;
    int *ggcslmem1,*ggcslmem2,*ggcslmem3,*ggcslmem4,*ggcslmem4a;
    int ggcslcount1,ggcslcount2,ggcslcount3,ggcslcount4,ggcslcount4a;
    int ggcslsize1,ggcslsize2,ggcslsize3,ggcslsize4,ggcslsize4a;

    
    // SLICE parallel
	int** gcslpara1;
	int** gcslpara2;
	int** gcslpara3;
	int** gcslpara4;

	int** gcslparaco1;
	int** gcslparaco2;
	int** gcslparaco3;
	int** gcslparaco4;
    
    
    // flow parameters
    const double cmu;
    double deltax,sigT,Ui,Ua,Uo;

    // Hydrodynamics Models
    int A10;
    
    // SFLOW
	int A209,A210,A211,A212,A214,A215,A216,A217,A218,A219,A220,A221,A230,A240,A241,A242,A243,A244,A245,A246,A248;
    int A251,A260;
    double A261,A262;
    double A223,A244_val,A245_val,A247,A249,A251_val;
    double A250;
    
    // FNPF
    int A310,A311,A312,A313,A320,A321,A322,A323,A329,A343,A344,A345,A347,A348;
    double A340,A341,A342,A344_val,A345_val,A346;
    int A350,A351,A352,A353,A357,A358,A361,A362,A363,A368;
    double A354,A355,A356,A365; 
    
    // NSEWAVE
    int A410;
    double A440;
    
    // NHFLOW
    int A501,A510,A511,A512,A514,A515,A516,A517,A518;
    int A520,A521;
    double A522,A523;
    double A531;
    int A540,A543;
    double A541,A542,A544,A545;
    int A550,A551,A552,A553;
    int A560;
    double A560_xs,A560_xe,A560_ys,A560_ye;
    int A561;
    double *A561_xs,*A561_xe,*A561_ys,*A561_ye,*A561_zs,*A561_ze;
    int A564;
    double *A564_xc,*A564_yc,*A564_zs,*A564_ze,*A564_r;
    
	// boundary conditions
	int B10,B20,B23;
    int B30,B32,B33;
    double B31,B32_x,B32_y,B32_z;    
    int B60,B61,B71,B75,B76,B77,B84,B85,B81,B82,B86,B87,B89,B90,B91,B92,B93,B94,B98,B99,B101,B105,B106,B107;
	int B136,B138,B138_1,B138_2,B139;
    int B180,B191,B192,B240,B241,B242,B243;
	double B29,B50,B51,B52,B53,B54,B55,B56,B81_1,B81_2,B81_3,B83,B117,B87_1,B87_2,B88;
	double B91_1,B91_2,B93_1,B93_2,B94_wdt,B96_1,B96_2,B102,B105_1,B105_2,B105_3;
	double *B71_val,*B71_dist,*B71_b,*B71_x,*B71_y;
	double *B106_b,*B106_x,*B106_y;
    double *B107_xs,*B107_xe,*B107_ys, *B107_ye, *B107_d;
    int B108;
    double *B108_xs,*B108_xe,*B108_ys, *B108_ye, *B108_d;
    int B110;
    double B110_zs,B110_ze;
	double B111_zs,B111_ze;
    double B112_zs,B112_z2,B112_ze;
    int B115,B116,B125,B127;
    double B120,B122,B123,B125_y;
    int B130,B133;
    double B131,B132_s,B132_e;
    double B134,B135;
    int B160, B170;
    int B181,B182,B183;
	double B181_1,B181_2,B181_3,B182_1,B182_2,B182_3,B183_1,B183_2,B183_3;
	double B191_1,B191_2,B191_3,B191_4,B192_1,B192_2,B192_3,B192_4;
	double B194_s,B194_e;
    
    int B411,B412,B413,B414,B415,B416,B417,B418,B421,B422;
    int *B411_ID;
    double *B411_Q;
    int *B412_ID;
    double *B412_pressBC;
    int *B413_ID;
    double *B413_h;
    int *B414_ID;
    double *B414_Uio;
    int *B415_ID;
    double *B415_U,*B415_V,*B415_W;
    int *B416_ID;
    double *B416_alpha;
    int *B417_ID;
    double *B417_Nx,*B417_Ny,*B417_Nz;
    int *B418_ID;
    int *B418_pio;
    int *B421_ID,*B421_Q;
    int *B422_ID,*B422_FSF;
    int B440;
    int *B440_ID,*B440_face;
    double *B440_xs,*B440_xe,*B440_ys,*B440_ye;
    int B441;
    int *B441_ID,*B441_face;
    double *B441_xs,*B441_xe,*B441_ys,*B441_ye,*B441_zs,*B441_ze;
    int B442;
    int *B442_ID,*B442_face;
    double *B442_xm,*B442_ym,*B442_zm,*B442_r;
    
	double *B240_D, *B240_C, *B240_xs, *B240_xe, *B240_ys, *B240_ye, *B240_zs, *B240_ze;
    double B260,B264,B267;
    int B269,B270;
    double *B270_xs, *B270_xe, *B270_ys, *B270_ye, *B270_zs, *B270_ze, *B270_n, *B270_d50, *B270_alpha, *B270_beta;
    int B274;
    double *B274_xc,*B274_yc,*B274_zs,*B274_ze,*B274_r, *B274_n, *B274_d50, *B274_alpha, *B274_beta;
    int B281;
    double *B281_xs, *B281_xe, *B281_ys, *B281_ye, *B281_zs, *B281_ze, *B281_n, *B281_d50, *B281_alpha, *B281_beta;
    int B282;
    double *B282_xs, *B282_xe, *B282_ys, *B282_ye, *B282_zs, *B282_ze, *B282_n, *B282_d50, *B282_alpha, *B282_beta;
	int B291;
    double *B291_xs, *B291_xe, *B291_ys, *B291_ye, *B291_zs, *B291_ze, *B291_d, *B291_n, *B291_d50, *B291_alpha, *B291_beta;
    int B295;
    int B308,B310,B311;
    double B309;
    double *B310_xs, *B310_xe, *B310_ys, *B310_ye, *B310_zs, *B310_ze, *B310_N, *B310_D, *B310_Cd;
    double *B311_xm, *B311_ym, *B311_r, *B311_zs, *B311ze, *B311_N, *B311_D, *B311_Cd;
    int B321;
    double *B321_xs, *B321_xe, *B321_ys, *B321_ye, *B321_zs, *B321_ze, *B321_N, *B321_D, *B321_Cd;
    int B322;
    double *B322_xs, *B322_xe, *B322_ys, *B322_ye, *B322_zs, *B322_ze, *B322_N, *B322_D, *B322_Cd;
	
    // Concentration Options
    double C1; //density concentration in water
    double C2; //viscosity water + concentration
    double C3; //density concentration in air
    double C4; //viscosity air + concentration
    double C5; //Schmidt number
    int C9; //only phase 1 concentratio
    int C10; //concentration transfer on/off
    int C15; //concentration convectio
    int C20; //concentration diffusio
    double C50_1; //fill ration concentration area 1
    double C50_2; //fill ration concentration area 2
    double C51; //i-dir zero level set start
    double C52; //j-dir zero level set start
    double C53; //k-dir zero level set start
    double C54; //i-dir zero level set en
    double C55; //j-dir zero level set en
    double C56; //k-dir zero level set en
    double C57_1; //a, plan
    double C57_2; //b
    double C57_3; //c
    double C57_4; //d
    double C58_1; //x0, spher
    double C58_2; //y0
    double C58_3; //z0
    double C58_4; //r
    int C75; //number of tiltboxes
    double *C75_x;
    double *C75_z;
    double *C75_a;
    double *C75_s;
    double *C75_l;
    double *C75_v;

	// discretization
    int D10; //convection scheme
    int D11; //convection velocity scheme
    int D20; //diffusion scheme
    int D21; //print out implicit diffusion time and iterations
    int D30; //pressure scheme
    int D31; //normalize pressure to free surface
    int D33; //corner cells sigma grid Poisson matrix
    int D37; //type of FSFBC for single fluid flow

	// Free Surface
    int F10; //free surface scheme
    int F30; //level set scheme
    int F31; //particle level se
    int F32; //number of particles per cell
    int F34; //printout iteration for pls
    int F35; //convection scheme for fsf
    int F36; //RK3 scheme
    int F40; //reini scheme
    int F44; //number reini time step
    int F46; //picard iteration for lsm or re
    int F47; //number of picard iterations
    int F49; //no reinitialization for interface nodes
    int F50; //bc phi, 1: inflow or 2: outflow
    int F150; //benchmark
    int F151; //benchmark inverse sign of level se
    double F33; //factor for pls vec allocation
    double F39; //reini constraint relaxation factor
    double F42; //maxlength
    double F43; //factor for reini timestep
    double F45; //factor for calculation of epsi
    double F51; //i-dir zero level set start
    double F52; //j-dir zero level set start
    double F53; //k-dir zero level set start
    double F54; //i-dir zero level set en
    double F55; //j-dir zero level set en
    double F56; //k-dir zero level set en
    int F50_flag; //flag for lsm descriptio
    double F57_1; //a, plan
    double F57_2; //b
    double F57_3; //c
    double F57_4; //d
    double F58_1; //x0, spher
    double F58_2; //y0
    double F58_3; //z0
    double F58_4; //r
    double F59_xm; //xm, cylinder
    double F59_ym; //ym
    double F59_zs; //zs
    double F59_ze; //z
    double F59_r; //r
    double F60; //ini z-dir
    double F61; //inflow  ini
    double F62; //outflow  ini
    double F63; //xstart phi interpolate with outflow h
    int F64; //fsf plane with angle on/off
    double F64_xs; //xs
    double F64_ys; //xs
    double F64_zs; //xs
    double F64_alpha; //alpha
    int F70; //number of phi 1 ini boxes
    double *F70_xs;
    double *F70_xe;
    double *F70_ys;
    double *F70_ye;
    double *F70_zs;
    double *F70_ze;
    int F71; //number of phi 2 ini boxes
    double *F71_xs;
    double *F71_xe;
    double *F71_ys;
    double *F71_ye;
    double *F71_zs;
    double *F71_ze;
    int F72; //number of phi 1 ini regions
    double *F72_xs;
    double *F72_xe;
    double *F72_ys;
    double *F72_ye;
    double *F72_h;
    int F80; //time scheme VOF
    int F85; //convection scheme VOF
    double F84; //cgamma for vof compression
    int F300; //multiphase flow level se
    int F305; //multiphase flow lsm convectio
    int F310; //multiphase flow re
    int F350; //multiphase flow fix level set inflow/outflow
    double F321; //epsi12
    double F322; //epsi13
    double F323; //epsi23
    double F360; //ini x-dir ls1
    double F361; //ini y-dir ls1
    double F362; //ini z-dir ls1
    int F369; //number of phi 1 ini tiltboxes ls1
    int F370; //number of phi 1 ini boxes ls1
    int F371; //number of phi 2 ini boxes ls1
    int F374; //number of pos ls1 ycyl
    int F375; //number of neg ls1 ycyl
    int F378; //number of pos ls1 sphere
    int F379; //number of neg ls1 sphere
    double *F369_x;
    double *F369_z;
    double *F369_a;
    double *F369_s;
    double *F369_l;
    double *F369_v;
    double *F370_xs;
    double *F370_xe;
    double *F370_ys;
    double *F370_ye;
    double *F370_zs;
    double *F370_ze;
    double *F371_xs;
    double *F371_xe;
    double *F371_ys;
    double *F371_ye;
    double *F371_zs;
    double *F371_ze;
    double *F374_xc;
    double *F374_zc;
    double *F374_r;
    double *F375_xc;
    double *F375_zc;
    double *F375_r;
    double *F378_xc;
    double *F378_yc;
    double *F378_zc;
    double *F378_r;
    double *F379_xc;
    double *F379_yc;
    double *F379_zc;
    double *F379_r;
    double F380; //ini x-dir ls2
    double F381; //ini y-dir ls2
    double F382; //ini z-dir ls2
    int F390; //number of phi 1 ini boxes ls2
    int F391; //number of phi 2 ini boxes ls2
    int F394; //number of pos ls2 ycyl
    int F395; //number of neg ls2 ycyl
    int F398; //number of pos ls2 sphere
    int F399; //number of neg ls2 sphere
    double *F390_xs;
    double *F390_xe;
    double *F390_ys;
    double *F390_ye;
    double *F390_zs;
    double *F390_ze;
    double *F391_xs;
    double *F391_xe;
    double *F391_ys;
    double *F391_ye;
    double *F391_zs;
    double *F391_ze;
    double *F394_xc;
    double *F394_zc;
    double *F394_r;
    double *F395_xc;
    double *F395_zc;
    double *F395_r;
    double *F398_xc;
    double *F398_yc;
    double *F398_zc;
    double *F398_r;
    double *F399_xc;
    double *F399_yc;
    double *F399_zc;
    double *F399_r;
    
	// Grid Options
    int G1; //xmargin inflow
    int G2; //sigma grid
    int G3; //solid forcing
    int G10; //xmargin inflow
    int G11; //ymargin righ
    int G12; //zmargin bottom
    int G20; //xmargin outflow
    int G21; //ymargin lef
    int G22; //zmargin top
    int G30; //extrapolated ghost cells
    int G40; //reini scheme for topo

	// Heat Options
    double H1; //thermal diffusivity water
    double H2; //thermal diffusivity air
    int H3; //ype of density calculatio
    int H4; //use beta coeff
    int H9; //air-water assignme
    int H10; //heat transfer on/off
    int H15; //convection for heat transfer
    double H4_beta1; //beta1
    double H4_beta2; //beta2
    double H50_1; //temperature 1
    double H50_2; //temperature 2
    double H51; //i-dir zero level set start
    double H52; //j-dir zero level set start
    double H53; //k-dir zero level set start
    double H54; //i-dir zero level set en
    double H55; //j-dir zero level set en
    double H56; //k-dir zero level set en
    double H57_1; //a, plan
    double H57_2; //
    double H57_3; //c
    double H57_4; //
    double H58_1; //x0, spher
    double H58_2; //y0
    double H58_3; //z0
    double H58_4; //r
    int H61; //heat bc
    int H62; //heat bc
    int H63; //heat bc
    int H64; //heat bc
    int H65; //heat bc
    int H66; //heat bc
    double H61_T; //heat bc
    double H62_T; //heat bc
    double H63_T; //heat bc
    double H64_T; //heat bc
    double H65_T; //heat bc
    double H66_T; //heat bc
	
	// Initialize Options
    int I10; //initialize all
    int I11; //initialize velocities with potential flow
    int I12; //initialize pressure
    int I13; //initialize turbulence
    int I30; //Fully intialize NWT
    int I40; //ini from state file
    int I41; //ID of state file
    int I44; //FNPF state with F
    int I56; //pressure above F56 set to zero
    double I21; //int set phase 2 velocities to zero after potential flow solver
    double I50; //simtime ini
    double I55; //reference pressur
    double I58_1; //vertical velocity for sphere initialization
    double I58_2; //radius for sphere initialization
    int I230; //read 2D flowfile
    double I231; //starting x for flowfi
    double I232; //starting y for flowfi
    double I233; //starting z for flowfi
    int I240; //read flowfile
    double I241; //delta t for flowfi

	// Numerical Options
    int N10; //linear poisson solver
    int N11; //precondioner
    int N40; //time scheme
    int N45; //max outer iter
    int N46; //max number of solver iterations
    int N48; //adaptive timestepping
    int N60; //maximum iteration of pjm correctio
    double N41; //total tim
    double N43; //stopping criteria convection-diffusion
    double N44; //stopping criteria pressur
    double N47; //relaxation factor for time stepping
    double N49; //max timestep or fixed timesteps
    double N50; //int adaptive timestepping meth
    double N61; //stopping criteria velocities

	// MPI Options
	int M10; //number of MPI processes

	// Print options
    int P10; //print file type
    int P11; //log print frequency
    int P12; //terminal print frequency
    int P15; //print file numbering
    int P16; //add timestamp to paraview files
    int P18; //option for phi print ou
    int P20; //h iteration file printed
    int P21; //time averaged vtk file print ou
    int P23; //print test to vtk file
    int P24; //print density to vtk file
    int P25; //print solid to vtk file
    int P26; //print cbed and conc to vtk file
    int P27; //print topo to vtk file
    int P28; //print fb to vtk file
    int P29; //print walldist to vtk file
    int P35; //print for interval
    int P40; //print state file
    int P41; //print state file each ith iteratio
    int P43; //state print out selected area
    int P44; //print out 3D potential for FNPF
    int P45; //print into single or continous state file
    int P50; //wave theory wave gages
    int P51; //print out wsf
    int P52; //print out wsfline in x-dir
    int P53; //print out wsfline for wave theory
    int P54; //ith iteration wsfline file  print ou
    int P56; //print out wsf line in y-dir
    int P57; //add aditional info to WSF gage in FNPF
    int P58; //print wave time series
    int P59; //print breaking wave log FNPF
    int P61; //print point probes
    int P62; //print line probes
    int P63; //print depth averaged point probe
    int P64; //print pressure probes
    int P65; //print velocity probes
    int P66; //print velocity probes from wave theory
    int P71; //print viscosity to vtk file
    int P72; //print vof functio
    int P73; //print hx and hy for sflow vtp
    int P74; //unused
    int P75; //print out vorticity vec
    int P76; //print out bedload
    int P77; //print out sediment parameters: 1
    int P78; //print out sediment parameters: 2
    int P79; //print out bed shear stress when running sediment transpor
    int P81; //force print ou
    int P82; //add eddyv to viscous force
    int P85; //ALE force print out for FNPF
    int P92; //force from water or from water+air
    int P101; //print sloshing forces
    int P120; //sediment log print ou
    int P121; //bed level gages
    int P122; //max bed level gages
    int P123; //topoline in x-directio
    int P124; //topoline in y-directio
    int P125; //bed shear stress gages
    int P126; //bed shear stress maxval
    int P140; //runup gage cylinder
    int P150; //number of data points to read from grid file
    int P151; //type of data
    int P152; //type of boundary condition for data
    int P166; //print discharge to terminal
    int P167; //discharge gages in x-directio
    int P168; //discharge gages in x-directio
    int P180; //print fsf
    int P181; //ith iteration fsf printed
    int P184; //time between file printout in iterations
    int P185; //time between file printout in seconds
    int P190; //print topo
    int P191; //ith iteration topo printed
    int P194; //time between file printout in iterations
    int P195; //time between file printout in seconds
    int P351; //print out wsf lsm1
    int P352; //print out wsf lsm2
    double P22; //start averging after transients
    double P30; //time between file printout in seconds
    double P34; //time between file printout in seconds for sediment
    double P42; //print state file each ith sec
    double *P35_ts;
    double *P35_te;
    double *P35_dt;
    double P43_xs;
    double P43_xe;
    double P43_ys;
    double P43_ye;
    int P46; //print state iteration window
    int P46_is;
    int P46_ie;
    int P47; //print state time window
    int P47_ts;
    int P47_te;
    double *P50_x;
    double *P50_y;
    double *P51_x;
    double *P51_y;
    double *P52_y;
    double *P56_x;
    double P55; //ith second wsfline files print out
    double *P58_x;
    double *P58_y;
    double *P58_T;
    double *P61_x;
    double *P61_y;
    double *P61_z;
    double *P62_xs;
    double *P62_ys;
    double *P62_zs;
    double *P62_xe;
    double *P62_ye;
    double *P62_ze;
    double *P63_x;
    double *P63_y;
    double *P64_x;
    double *P64_y;
    double *P64_z;
    double *P65_x;
    double *P65_y;
    double *P65_z;
    double *P66_x;
    double *P66_y;
    double *P66_z;
    double *P81_xs;
    double *P81_xe;
    double *P81_ys;
    double *P81_ye;
    double *P81_zs;
    double *P81_ze;
    double *P85_x;
    double *P85_y;
    double *P85_r;
    double *P85_cd;
    double *P85_cm;
    double P91; //factor used in force calculation algorithm
    double P101_xm;
    double P101_ym;
    double P101_zs;
    double P101_ze;
    double P101_r1;
    double P101_r2;
    int P110; //print significant wave heigh
    double P111; //start averging after transients
    double *P121_x;
    double *P121_y;
    double *P123_y;
    double *P124_x;
    double *P125_x;
    double *P125_y;
    int P131; //max wetdry in vtp
    int P132; //max wetdry as file
    int P133; //runup gage x-crossectio
    int P134; //runup gage y-crossectio
    double *P133_y;
    double *P134_y;
    double *P140_x;
    double *P140_y;
    double *P167_x;
    double P141; //int runup cylinder radius
    double *P168_x;
    double *P168_zs;
    double *P168_ze;
    double P182; //time between fsf file printout in seconds
    int *P184_its;
    int *P184_ite;
    int *P184_dit;
    double *P185_ts;
    double *P185_te;
    double *P185_dt;
    double P192; //time between topo file printout in seconds
    int *P194_its;
    int *P194_ite;
    int *P194_dit;
    double *P195_ts;
    double *P195_te;
    double *P195_dt;
    int P230; //print flowfile
    int P240; //print potentialfile
    double *P230_x;
    double *P240_x;
    double *P351_x;
    double *P351_y;
    double *P352_x;
    double *P352_y;
    
    // Particles
    int Q10,Q24,Q29,Q43;
    double Q21,Q22,Q23,Q25;
    double Q31;
    double Q41;
    
    int Q101,Q110;
    double *Q110_xs,*Q110_xe,*Q110_ys,*Q110_ye,*Q110_zs,*Q110_ze;
    int Q111,Q112,Q113;
    double Q111_x,Q112_y,Q113_z;
    
    int Q180,Q181;
    double Q182;
    

	// Sediment Transport
    int S10; //sediment transport module
    int S11; //bedload formula
    int S12; //Suspended Sediment, formula for boundary condition
    int S15; //synchronize sediment time step with main solver
    int S16; //bed shear stress formulatio
    int S17; //non-equillibrium bedload 
    int S25; //automatic sediment fall velocity
    int S27; //number of inner iterations
    int S32; //exner discretizatio
    int S33; //type of near bead velocity interpolatio
    int S34; //type of suspedned load D and E calculatio
    int S37; //number reini time step
    int S41; //type of sediment start criterio
    int S42; //type of sediment interval criterio
    int S43; //number of water iteration, before sediment transport starts
    int S44; //number of water timesteps between bed calculatio
    int S50; //bc phi, 1: inflow fix or 2: outflow fix, 3: both fix
    int S60; //time stepping for suspended sediments
    int S73; //distance for use relaxation method for the sediment bed
    int S77; //active sediment domain in x-directio
    int S78; //inflow guard
    int S79; //outflow guard
    int S80; //type of slope reductio
    int S83; //type of bedslope calc
    int S84; //type of critical bed shear stress reduction limiters
    int S90; //sandslide on/off
    int S91; //number of sandslide iterations
    int S100; //number of bed filter outer iterations
    int S101; //number of bed filter inner iterations
    double S13; //timestep for sediment transport
    double S14; //relaxation timestep size for sediment transport
    double S19; //total time sediment
    double S20; //sediment d50
    double S21; //factor for d50 for calculation of ks in bedshear routin
    double S22; //sediment density
    double S23; //sediment fall velocity
    double S24; //porosity of sediment layer
    double S26_a; //alpha for VRANS sediment
    double S26_b; //beta for VRANS sediment
    double S30; //Shields parameter
    double S45; //flow simulation time, before sediment transport starts
    double S46; //flow simulation time between bed calculation
    double S47; //t/T, before sediment transport starts
    double S48; //int nt/T between bed calculation
    double S57; //ini z-dir
    double S71; //int x start of erosion
    double S72; //int x end of erosion
    double S81; //midphi for slope reduction
    double S82; //delta phi for slope reduction
    double S92; //sandslide correction factor
    double S93; //delta phi for sandlide correciton
    double *S73_val;
    double *S73_dist;
    double *S73_b;
    double *S73_x;
    double *S73_y;
    double S77_xs; //active sediment domain x_start
    double S77_xe; //active sediment domain x_en

	// Turbulence
    int T10; //turbulence model
    int T12; //convection scheme
    int T21; //type of LES filter
    int T33; //kin source
    int T36; //explciti free surface dampong through dissipatio
    int T39; //blend fsf eddyv with sgs-eddyv
    int T41; //RANS stabilizatio
    int T44; //buouncy term
    double T31; //factor for limiter for eddy limiter in phase 1
    double T32; //factor for limiter for eddy limiter in phase 2
    double T35; //factor for limiter for eddy limiter near wa
    double T37; //int damping coefficient for T36
    double T38; //epsi fsf turbulence damping
    double T42; //lambda1 factor
    double T43; //komega wall BC velocity factor

	// Waterflow
    double W1; //density water
    double W2; //viscosity water
    double W3; //density air
    double W4; //viscosity air
    double W5; //surface tension between phase 1 and phase 2
    double W6; //density oi
    double W7; //viscosity oi
    double W10; //discharg
    double W_fb; //density of floating body
    int W11; //velocity inlet face 1
    int W12; //velocity inlet face 2
    int W13; //velocity inlet face 3
    int W14; //velocity inlet face 4
    int W15; //velocity inlet face 5
    int W16; //velocity inlet face 6
    double W11_u; //u-velocity inlet face 1
    double W11_v; //v-velocity inlet face 1
    double W11_w; //w-velocity inlet face 1
    double W12_u; //u-velocity inlet face 2
    double W12_v; //v-velocity inlet face 2
    double W12_w; //w-velocity inlet face 2
    double W13_u; //u-velocity inlet face 3
    double W13_v; //v-velocity inlet face 3
    double W13_w; //w-velocity inlet face 3
    double W14_u; //u-velocity inlet face 4
    double W14_v; //v-velocity inlet face 4
    double W14_w; //w-velocity inlet face 4
    double W15_u; //u-velocity inlet face 5
    double W15_v; //v-velocity inlet face 5
    double W15_w; //w-velocity inlet face 5
    double W16_u; //u-velocity inlet face 6
    double W16_v; //v-velocity inlet face 6
    double W16_w; //w-velocity inlet face 6
    double W20; //gi
    double W21; //gj
    double W22; //gk
    double W31; //temperature for air compressibility in celsius
    double W29_x; //pressure gradient x-direction
    double W29_y; //pressure gradient y-direction
    double W29_z; //pressure gradient z-direction
    int W30; //air compressibility on/off
    int W41; //velocity source phase 1
    double *W41_xc;
    double *W41_yc;
    double *W41_zs;
    double *W41_ze;
    double *W41_vel;
    double *W41_beta;
    double W50; //air inflow
    int W50_air; //air inflow switch
    int W90; //non-newtownian flow
    double W95; //nu_0
    double W96; //tau_0
    double W97; //K
    double W98; //int
    int W101; //turn on Mohr-Coloumb
    double W102_c; //c factor
    double W102_phi; //angle of repos
    double W103; //MC transition factor
    double W104; //shear rate dependent excess pore pressure factor
    int W110; //add rheology as source term or viscosity
    int W111; //which pressure for MC
    double W112; //threshold factor for pressure blening in W111 3
    
    // 6DOF
	double ufb,vfb,wfb;
	double pfb,qfb,rfb;
	double ufbi,vfbi,wfbi;
	double pfbi,qfbi,rfbi;
	double xg,yg,zg;
	double xgn,ygn,zgn;
	double phi_fb,theta_fb,psi_fb;
	double ufbmax, vfbmax, wfbmax;
    int X10,X12,X14,X15,X19,X11_u,X11_v,X11_w,X11_p,X11_q,X11_r,X21,X22,X23,X24,X31,X32,X33,X34,X38;
    int X39,X40,X45,X46,X47,X48,X49,X50,X60,X110,X120,X131,X132,X133;
	int X100,X101,X102,X103,X141,X142,X143,X153,X180,X181,X182,X183,X210,X211;
	int X310, X311, X312, X313, X314, X315, X320, X321, mooring_count, net_count;
	double X21_d,X22_m;
	double X23_x,X23_y,X23_z;
	double X24_Ix,X24_Iy,X24_Iz;	
	double X25_Cp,X25_Cq,X25_Cr;	
    double X26_Cu,X26_Cv,X26_Cw;	
	double X41,X42,X43,X44;
	double X100_x,X100_y,X100_z;
	double X101_phi, X101_theta, X101_psi;
	double X102_u, X102_v, X102_w;
	double X103_p, X103_q, X103_r;
	double *X110_xs,*X110_xe,*X110_ys,*X110_ye,*X110_zs,*X110_ze;
	double X120_rad,X120_xc,X120_yc,X120_zc;
	double X131_rad,X131_h,X131_xc,X131_yc,X131_zc;
	double X132_rad,X132_h,X132_xc,X132_yc,X132_zc;
	double X133_rad,X133_h,X133_xc,X133_yc,X133_zc;
	double X153_xs,X153_xe,X153_ys,X153_ye,X153_zs,X153_ze;
    int X163;
    double *X163_x1,*X163_y1,*X163_z1;
    double *X163_x2,*X163_y2,*X163_z2;
    double *X163_x3,*X163_y3,*X163_z3;
    double *X163_x4,*X163_y4,*X163_z4;
    double *X163_x5,*X163_y5,*X163_z5;
    double *X163_x6,*X163_y6,*X163_z6;
    int X164;
    double *X164_x1,*X164_y1,*X164_z1;
    double *X164_x2,*X164_y2,*X164_z2;
    double *X164_x3,*X164_y3,*X164_z3;
    double *X164_x4,*X164_y4,*X164_z4;
    double *X164_x5,*X164_y5,*X164_z5;
    double *X164_x6,*X164_y6,*X164_z6;
    double *X164_x7,*X164_y7,*X164_z7;
    double *X164_x8,*X164_y8,*X164_z8;
    double X181_x,X181_y,X181_z;
    double X182_x,X182_y,X182_z;
    double X183_x,X183_y,X183_z,X183_phi,X183_theta,X183_psi;
    int X185,X188;
    double X186;
    
    int X205;
    int X206,X207;
    double X206_ts,X206_te,X207_ts,X207_te;
	double X210_u,X210_v,X210_w;
	double X211_p,X211_q,X211_r;
    int X240;
    double X241,X242_x,X242_y,X242_z,X243;
    double *X311_xs,*X311_xe,*X311_ys,*X311_ye,*X311_zs,*X311_ze;
    double *X311_w,*X311_rho_c,*X311_EA,*X311_d,*X311_l,*X311_H,*X311_P,*X311_facT;
    double *X312_k,*X312_T0;
    double *X314_T, *X315_t;
    int *X320_type;
	double *X321_Sn,*X321_d,*X321_lambda,*X321_dk,*X321_rho,*X321_nd,*X321_nl;
    double *X322_D,*X322_L,*X322_x0,*X322_y0,*X322_z0,*X322_phi,*X322_theta,*X322_psi;
    int X324;
    double X323_m,X323_d,X323_l;
    double *X324_x,*X324_y,*X324_z;
    double X325_dt,X325_relX,X325_relY,X325_relZ;
    int X400;
    double X401_p0,X401_cl,X401_cb,X401_a;

    // FSI
    int Z10,Z11,FSI_count;
    double *Z11_x,*Z11_y,*Z11_z,*Z11_l,*Z11_w,*Z11_t,*Z11_rho,*Z11_e,*Z11_ix,*Z11_iy,*Z11_iz,*Z11_nu,*Z11_n;
    double Z12_ckx,Z12_cky,Z12_ckz,Z12_cdx,Z12_cdy,Z12_cdz;
	
	// Grid
	int Y40,Y50,Y60,Y71,Y72,Y73,Y74;

    // Test options
    int Y1,Y2,Y3,Y4,Y5;

	// time + iterations
	int inneriter,count,solveriter,preconiter,count_statestart;
    int solver_status;
    int sediter;
    double final_res;
	double dt,dt_old,simtime,viscmax;
	double mindt,maxdt;
	double umax,vmax,wmax,epsmax,kinmax,pressmin,pressmax,omegamax;
	double presstime,veltime,reinitime,turbtime,plstime,itertime;
	double sedsimtime,sedwavetime;
	double wavetime;
	double meantime,totaltime;
	double gcmeantime,gctotaltime;
	double Xmeantime,Xtotaltime;
	double maxbed, minbed;
	double susptime,maxtopovel;
	double gctime, xtime;
	double volume1,volume2,volume3;
	double Qi,Qo;
	double dtsed,sedtime,slidecells;
	double bedmax,bedmin;
	double field4time;
    double printtime, sedprinttime,fsfprinttime,probeprinttime,stateprinttime,exportprinttime;
    double partprinttime;

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
	double kintime,epstime;
	double poissontime, laplacetime;
    double sftime,fbtime,fsitime;
    double fbdt,fbmax;
    double sfdt,sfmax;
	double lsmtime,heattime,concentrationtime;
	double printouttime;
	double phimean,phiout,phiin;
    double fsfin,fsfout;
    double fsfinval,fsfoutval;
	double pcnorm,ucnorm,vcnorm,wcnorm;
    double alpha;
    double pressgage;

	// maxcoor
	double xcoormax,xcoormin,ycoormax,ycoormin,zcoormax,zcoormin;
	double maxlength;


	// wave coefficients
	double wT,wV,wH,wA,wL,wd,ww,wk,wC;
	double wHs,wAs,wwp,ww_s,ww_e,wTp;
	int wN;
    double wts,wte;
    
    // free surface
    double psi,psi0;
	int pressval;

// Boundary
    //int **boundary;
    int **fgc;

	static int knox,knoy,knoz;
	static int margin;

	static int xtp,ytp,ztp;
	static int xmax,ymax,zmax;

// PARALELL
    int mpirank;
	int gcx_1range1[7],gcx_3range1[7];
	int gcx_1range2[7],gcx_3range2[7];
	int gcx_1range3[7],gcx_3range3[7];
	int gcx_1range4[7],gcx_3range4[7];
	
// Non-Uniform Mesh    
    double *XN,*YN,*ZN;
    double *XP,*YP,*ZP;
    double *DXN,*DYN,*DZN;
    double *DXP,*DYP,*DZP;
    double *ZSN,*ZSP;
    double DXM,DYD,DXD;
    double DYM,DZM;
    
    double *RN,*SN,*TN;
    double *RP,*SP,*TP;
    double *DRN,*DSN,*DTN;
    double *DRP,*DSP,*DTP;
    double DRM,DSM,DTM;
    double DX,DY,DZ;
    double *DRDXN,*DSDYN,*DTDZN;
    double *DRDXP,*DSDYP,*DTDZP;
    double *DDRDDXN,*DDSDDYN,*DDTDDZN;
    double *DDRDDXP,*DDSDDYP,*DDTDDZP;
    
    weno_nug_func *wenofunc;
    
// sigma coordinate
    double *sig;
    double *sigx,*sigy,*sigz,*sigt;
    double *sigxx;
    
private:
	void clear(char&, int&);
    
};

#endif
