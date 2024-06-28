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
    int A10; //turn on wave models
    
    // SFLOW
    int A209; //interpolation sweeps for bed
    int A210; //time scheme for SFLOW velocities
    int A211; //convection scheme for SLOW velocities
    int A212; //diffusion treatment for SLOW velocities
    int A214; //convection for vertical velocity
    int A215; //conservative discretization
    int A216; //convection velocity
    int A217; //slip or no-slip boundary conditions
    int A218; //turn on roughness
    int A219; //additional courant number constra
    int A220; //non-hydrostatic pressure scheme for SFLOW
    int A221; //hydrostatic pressure scheme for SFLOW
    int A230; //turn on Boussinesq wave model
    int A240; //FSF algorithm SFLOW
    int A241; //discretization of water level SFLOW
    int A242; //hydostatic pressure for shallow areas
    int A243; //turn on wetting-drying
    int A244; //double absolute wetting criterion value
    int A245; //dx-based relative wetting citerion
    int A246; //turn on breaking
    int A248; //turn on breaking persistence
    int A251; //double fsf-slope in x-dir
    int A260; //turbulence model
    double A261; //length scale factor
    double A262; //parabolic turbulence model factor
    double A223; //blending factor hydrostatic pressure gradient
    double A244_val; //absolute wetting criterion va
    double A245_val; //dx-based relative wetting citerion va
    double A247; //breaking parameter alpha
    double A249; //breaking persistence parameter beta
    double A251_val;
    double A250; //viscosity breaking wav
    
    // FNPF
    int A310; //time scheme for FNPF velocities
    int A311; //convection scheme for FNPF velocities
    int A312; //discretization for second-order gradie
    int A313; //discretization for bed bc
    int A320; //order of Laplace equation
    int A321; //boundary condition order for 4th-order Laplace equation
    int A322; //maxiter for 4th-order Laplace after 2nd-order solution
    int A323; //PTF FSF extrapolation
    int A329; //wave maker BC order
    int A343; //turn on wetting-drying
    int A344; //absolute wetting criterion
    int A345; //dx-based relative wetting citerion
    int A347; //coastline relaxation for Fi and eta
    int A348; //beach relaxation for Fi and eta
    double A340; //minimum water depth
    double A341; //coastline damping distance factor for dxm
    double A342; //coastline damping absolute distanc
    double A344_val; //absolute wetting criterion va
    double A345_val; //dx-based relative wetting citerion va
    double A346; //viscosity damping within the coastlin
    int A350; //turn on breaking (which method)
    int A351; //type of breaking detection (deep / shallow)
    int A352; //additional filtering to viscosity based breaking
    int A353; //breaking wave identification algorithm
    int A357; //breaking for Fi and eta
    int A358; //breaking algorithm version
    int A361; //breaking filter outer iter
    int A362; //breaking filter inner iter
    int A363; //breaking filter width
    int A368; //breaking waves in numerical beach
    double A354; //breaking parameter alpha
    double A355; //breaking parameter slope alpha
    double A356; //breaking parameter slope beta
    double A365; //viscosity breaking wav
    
    // NSEWAVE
    int A410; //scheme eta
    double A440; //epsi for depth integration
    
    // NHFLOW
    int A501; //nhf mode
    int A510; //NFHLOW time scheme
    int A511; //NHFLOW HLL scheme
    int A512; //NHFLOW diffusion
    int A514; //NHFLOW reconstruction 
    int A515; //NHFLOW KFSFBC scheme
    int A516; //NFHLOW KFSFBED scheme
    int A517; //NHFLOW omega_sig scheme
    int A518; //NHFLOW bed BC
    int A520; //NFHLOW non-hydrostatic pressure scheme
    int A521; //unused
    double A522; //p_alpha
    double A523; //p_gamma
    double A531; //Froude number limiter
    int A540; //NFHLOW fsf scheme
    int A543; //NHFLOW wetting & drying or coastline
    double A541; //coastline damping distance factor for dxm
    double A542; //coastline damping absolute distanc
    double A544; //wetting & drying criterion
    double A545; //deep criterion
    int A550; //turn on breaking (which method)
    int A551; //type of breaking detection (deep / shallow)
    int A552; //additional filtering to viscosity based breaking
    int A553; //breaking in very shallow regions turned onf
    int A560; //block eta
    double A560_xs;
    double A560_xe;
    double A560_ys;
    double A560_ye;
    int A561; //solid box
    double *A561_xs;
    double *A561_xe;
    double *A561_ys;
    double *A561_ye;
    double *A561_zs;
    double *A561_ze;
    int A564; //solid vertical cylinder
    double *A564_xc;
    double *A564_yc;
    double *A564_zs;
    double *A564_ze;
    double *A564_r;
    
	// boundary conditions
    int B10; //wall laws velocities on/off
    int B20; //slip or no-slip boundary condition for velocity
    int B23; //ghostcell extrapolation or refective
    int B30; //type of pressure reference po
    int B32; //pressure reference location
    int B33; //pressure gage virtual or inline
    double B31; //pressure reference va
    double B32_x; //pressure reference location
    double B32_y; //pressure reference location
    double B32_z; //pressure reference location
    int B60; //ioflow discharge
    int B61; //plain or logarithmic inflow profile
    int B71; //double distance for use relaxation method for fixed water level 
    int B75; //type of outflow boundary conditions
    int B76; //type of pressure inlet boundary condition
    int B77; //outflow pressure controlled or free stream
    int B84; //Peak enhance method
    int B85; //PM or JONSWAP spectrum for irregular waves
    int B81; //focussed wave parameter
    int B82; //type of focus point and time calculation
    int B86; //number of regular waves for irregular wave generation
    int B87; //give ws and we for irregular wave generation
    int B89; //wave generation optimization
    int B90; //iowave
    int B91; //wave parameter wL
    int B92; //wave type
    int B93; //wave parameter wT
    int B94; //set water depth for wave theory
    int B98; //type of wave generation
    int B99; //type of numerical beach
    int B101; //ramp function wave geneartion
    int B105; //wave generation origin changed
    int B106; //read wave generation orig
    int B107; //read numerical beach orig
    int B136; //double summation method frequency vector
    int B138; //seed number multidir waves
    int B138_1;
    int B138_2;
    int B139; //seed number wave spectrum
    int B180; //gravity waves
    int B191; //rotation around x-axis
    int B192; //rotation around y-axis
    int B240; //porous media
    int B241; //porous media in x-direction
    int B242; //porous media in y-direction
    int B243; //porous media in z-direction
    double B29; //gamma for gc image point
    double B50; //global wall roughness ks
    double B51; //global wall roughness ks
    double B52; //global wall roughness ks
    double B53; //global wall roughness ks
    double B54; //global wall roughness ks
    double B55; //global wall roughness ks
    double B56; //global wall roughness ks
    double B81_1;
    double B81_2;
    double B81_3; //unidirectional focused wave y is condisered 0
    double B83; //wave steepness parameter for focused breaking waves
    double B117; //starting time shift for timeseries input
    double B87_1;
    double B87_2;
    double B88; //gamma for JONSWAP spectrum
    double B91_1; //wave amplit
    double B91_2; //wave length
    double B93_1; //wave amplit
    double B93_2; //wave peri
    double B94_wdt; //water depth for wave theory
    double B96_1; //dist1 for wave relax
    double B96_2; //dist2 for wave relax
    double B102; //factor ramp function wave generation
    double B105_1; //wave generation direction and line origin in Cartesian coordianate system
    double B105_2; //wave generation line x origin
    double B105_3; //wave generation line y origin
    double *B71_val;
    double *B71_dist;
    double *B71_b;
    double *B71_x;
    double *B71_y;
    double *B106_b;
    double *B106_x;
    double *B106_y;
    double *B107_xs;
    double *B107_xe;
    double *B107_ys;
    double *B107_ye;
    double *B107_d;
    int B108; //read wave generation  orig
    double *B108_xs;
    double *B108_xe;
    double *B108_ys;
    double *B108_ye;
    double *B108_d;
    int B110; //read wave generation  orig
    double B110_zs;
    double B110_ze;
    double B111_zs; //flap start
    double B111_ze; //flap en
    double B112_zs; //flap start
    double B112_z2; //flap2 end/flap2 start
    double B112_ze; //flag en
    int B115; //activate vertical velocity component for flap wavemaker theory
    int B116; //x or beta input for flap wavemaker theories
    int B125; //take 2D slice input for HDC
    int B127; //turn of y-dir velociteis for HDC 
    double B120; //delta t for wave generation
    double B122; //int air velocity on/off for active wave generation
    double B123; //flap AWA hinge location
    double B125_y; //2D slice y-coor input for HDC
    int B130; //directional spreading for irregular waves
    int B133; //number of direction intervals for spreading function
    double B131; //main direction for multidirectional irregular waves
    double B132_s; //start directional spreading
    double B132_e; //end directional spreading
    double B134; //shape parameter for spreading function
    double B135; //peak va
    int B160; //number of vertical layers for 2D wave generation
    int B170; //number of Fourier modes for the generation of steady surface gravity waves
    int B181; //x-dir motion
    int B182; //y-dir motion
    int B183; //z-dir motion
    double B181_1; //x-acceleration amplit
    double B181_2; //x-acceleration frequency
    double B181_3; //wave phase chang
    double B182_1; //y-acceleration amplit
    double B182_2; //y-acceleration frequency
    double B182_3; //wave phase chang
    double B183_1; //z-acceleration amplit
    double B183_2; //z-acceleration frequency
    double B183_3; //wave phase chang
    double B191_1; //angle for rotation around x-axis
    double B191_2; //frequency for rotation around x-axis
    double B191_3; //y-coordinate for rotation around x-axis
    double B191_4; //z-coordinate for rotation around x-axis
    double B192_1; //angle for rotation around y-axis
    double B192_2; //frequency forrotation around y-axis
    double B192_3; //x-coordinate for rotation around y-axis
    double B192_4; //z-coordinate for rotation around y-axis
    double B194_s; //start rotation
    double B194_e; //end rotation
    int B411; //patchBC discharge
    int B412; //patchBC pressure BC
    int B413; //patchBC waterlevel
    int B414; //patchBC perpendicular velocity
    int B415; //patchBC velocity components
    int B416; //patchBC horizontal inflow angle
    int B417; //patchBC inflow normals
    int B418; //patchBC outflow pressure condition
    int B421; //patchBC hydrograph discharge
    int B422; //patchBC hydrograph waterlevel
    int *B411_ID;
    double *B411_Q;
    int *B412_ID;
    double *B412_pressBC;
    int *B413_ID;
    double *B413_h;
    int *B414_ID;
    double *B414_Uio;
    int *B415_ID;
    double *B415_U;
    double *B415_V;
    double *B415_W;
    int *B416_ID;
    double *B416_alpha;
    int *B417_ID;
    double *B417_Nx;
    double *B417_Ny;
    double *B417_Nz;
    int *B418_ID;
    int *B418_pio;
    int *B421_ID;
    int *B421_Q;
    int *B422_ID;
    int *B422_FSF;
    int B440; //patch BC inflow line
    int *B440_ID;
    int *B440_face;
    double *B440_xs;
    double *B440_xe;
    double *B440_ys;
    double *B440_ye;
    int B441; //rectangular inflow patch BC
    int *B441_ID;
    int *B441_face;
    double *B441_xs;
    double *B441_xe;
    double *B441_ys;
    double *B441_ye;
    double *B441_zs;
    double *B441_ze;
    int B442; //circular inflow patch BC
    int *B442_ID;
    int *B442_face;
    double *B442_xm;
    double *B442_ym;
    double *B442_zm;
    double *B442_r;
    double *B240_D;
    double *B240_C;
    double *B240_xs;
    double *B240_xe;
    double *B240_ys;
    double *B240_ye;
    double *B240_zs;
    double *B240_ze;
    double B260; //C coefficient for VRANS
    double B264; //KC number for VRANS
    double B267; //d50 for VRANS
    int B269; //VRANS on/off -> assigned as 1 for VRANS Structure, 2 for Vegetation, 3 for Net interaction
    int B270; //VRANS porous media box
    double *B270_xs;
    double *B270_xe;
    double *B270_ys;
    double *B270_ye;
    double *B270_zs;
    double *B270_ze;
    double *B270_n;
    double *B270_d50;
    double *B270_alpha;
    double *B270_beta;
    int B274; //VRANS porous media vertical cylinder
    double *B274_xc;
    double *B274_yc;
    double *B274_zs;
    double *B274_ze;
    double *B274_r;
    double *B274_n;
    double *B274_d50;
    double *B274_alpha;
    double *B274_beta;
    int B281; //VRANS porous media wedge in x-direction
    double *B281_xs;
    double *B281_xe;
    double *B281_ys;
    double *B281_ye;
    double *B281_zs;
    double *B281_ze;
    double *B281_n;
    double *B281_d50;
    double *B281_alpha;
    double *B281_beta;
    int B282; //VRANS porous media wedge in y-direction
    double *B282_xs;
    double *B282_xe;
    double *B282_ys;
    double *B282_ye;
    double *B282_zs;
    double *B282_ze;
    double *B282_n;
    double *B282_d50;
    double *B282_alpha;
    double *B282_beta;
    int B291; //VRANS porous media plate in x-direction
    double *B291_xs;
    double *B291_xe;
    double *B291_ys;
    double *B291_ye;
    double *B291_zs;
    double *B291_ze;
    double *B291_d;
    double *B291_n;
    double *B291_d50;
    double *B291_alpha;
    double *B291_beta;
    int B295;
    int B308; //porosity effects on fluid acceleration for vegetation
    int B310; //VRANS vegetation box
    int B311;
    double B309; //Cm for vegetation
    double *B310_xs;
    double *B310_xe;
    double *B310_ys;
    double *B310_ye;
    double *B310_zs;
    double *B310_ze;
    double *B310_N;
    double *B310_D;
    double *B310_Cd;
    double *B311_xm;
    double *B311_ym;
    double *B311_r;
    double *B311_zs;
    double *B311ze;
    double *B311_N;
    double *B311_D;
    double *B311_Cd;
    int B321; //VRANS vegetation wedge in x-direction
    double *B321_xs;
    double *B321_xe;
    double *B321_ys;
    double *B321_ye;
    double *B321_zs;
    double *B321_ze;
    double *B321_N;
    double *B321_D;
    double *B321_Cd;
    int B322; //VRANS vegetation wedge in y-direction
    double *B322_xs;
    double *B322_xe;
    double *B322_ys;
    double *B322_ye;
    double *B322_zs;
    double *B322_ze;
    double *B322_N;
    double *B322_D;
    double *B322_Cd;
	
    // Concentration Options
    double C1; //density concentration in water
    double C2; //viscosity water + concentration
    double C3; //density concentration in air
    double C4; //viscosity air + concentration
    double C5; //Schmidt number
    int C9; //only phase 1 concentration
    int C10; //concentration transfer on/off
    int C15; //concentration convection
    int C20; //concentration diffusion
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
    int F50_flag; //flag for lsm description
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
    int F305; //multiphase flow lsm convection
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
    int H3; //ype of density calculation
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
    int N60; //maximum iteration of pjm correction
    double N41; //total time
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
    int P18; //option for phi print out
    int P20; //h iteration file printed
    int P21; //time averaged vtk file print out
    int P23; //print test to vtk file
    int P24; //print density to vtk file
    int P25; //print solid to vtk file
    int P26; //print cbed and conc to vtk file
    int P27; //print topo to vtk file
    int P28; //print fb to vtk file
    int P29; //print walldist to vtk file
    int P35; //print for interval
    int P40; //print state file
    int P41; //print state file each ith iteration
    int P43; //state print out selected area
    int P44; //print out 3D potential for FNPF
    int P45; //print into single or continous state file
    int P50; //wave theory wave gages
    int P51; //print out wsf
    int P52; //print out wsfline in x-dir
    int P53; //print out wsfline for wave theory
    int P54; //ith iteration wsfline file  print out
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
    int P72; //print vof function
    int P73; //print hx and hy for sflow vtp
    int P74; //unused
    int P75; //print out vorticity vec
    int P76; //print out bedload
    int P77; //print out sediment parameters: 1
    int P78; //print out sediment parameters: 2
    int P79; //print out bed shear stress when running sediment transport
    int P81; //force print out
    int P82; //add eddyv to viscous force
    int P85; //ALE force print out for FNPF
    int P92; //force from water or from water+air
    int P101; //print sloshing forces
    int P120; //sediment log print out
    int P121; //bed level gages
    int P122; //max bed level gages
    int P123; //topoline in x-direction
    int P124; //topoline in y-direction
    int P125; //bed shear stress gages
    int P126; //bed shear stress maxval
    int P140; //runup gage cylinder
    int P150; //number of data points to read from grid file
    int P151; //type of data
    int P152; //type of boundary condition for data
    int P166; //print discharge to terminal
    int P167; //discharge gages in x-direction
    int P168; //discharge gages in x-direction
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
    int P133; //runup gage x-crossection
    int P134; //runup gage y-crossection
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
    int S16; //bed shear stress formulation
    int S17; //non-equillibrium bedload 
    int S25; //automatic sediment fall velocity
    int S27; //number of inner iterations
    int S32; //exner discretization
    int S33; //type of near bead velocity interpolation
    int S34; //type of suspedned load D and E calculation
    int S37; //number reini time step
    int S41; //type of sediment start criterion
    int S42; //type of sediment interval criterion
    int S43; //number of water iteration, before sediment transport starts
    int S44; //number of water timesteps between bed calculation
    int S50; //bc phi, 1: inflow fix or 2: outflow fix, 3: both fix
    int S60; //time stepping for suspended sediments
    int S73; //distance for use relaxation method for the sediment bed
    int S77; //active sediment domain in x-direction
    int S78; //inflow guard
    int S79; //outflow guard
    int S80; //type of slope reduction
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
    int T36; //explciti free surface dampong through dissipation
    int T39; //blend fsf eddyv with sgs-eddyv
    int T41; //RANS stabilization
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
    double ufb;
    double vfb;
    double wfb;
    double pfb;
    double qfb;
    double rfb;
    double ufbi;
    double vfbi;
    double wfbi;
    double pfbi;
    double qfbi;
    double rfbi;
    double xg;
    double yg;
    double zg;
    double xgn;
    double ygn;
    double zgn;
    double phi_fb;
    double theta_fb;
    double psi_fb;
    double ufbmax;
    double vfbmax;
    double wfbmax;
    int X10; //turn 6DOF o
    int X12; //turn force calculation o
    int X14; //tangential velocity 
    int X15; //density treatment for direct forcing
    int X19; //print out interval 6DOF log files
    int X11_u; //turn on degrees of freedom
    int X11_v; //turn on degrees of freedom
    int X11_w; //turn on degrees of freedom
    int X11_p; //turn on degrees of freedom
    int X11_q; //turn on degrees of freedom
    int X11_r; //turn on degrees of freedom
    int X21; //presribe homogeneous density floating body
    int X22; //prescribe mass floating body
    int X23; //prescribe center of gravity
    int X24; //prescribe moments of inertia
    int X31; //boundary conditions for parallel velocity on floating body
    int X32; //boundary conditions for orthogonal velocity on floating body
    int X33; //boundary conditions for pressure on floating body
    int X34; //boundary treatment for new solid velocity cells
    int X38;
    int X39; //type of viscous force calculation
    int X40; //type of force calculation
    int X45; //type of lsm convection disc at fb
    int X46; //density smoothing inside fb
    int X47; //reini diffusion inside fb
    int X48;
    int X49;
    int X50; //type of print out format for 6DOF structure
    int X60; //type of print of force calculation
    int X110; //rectangular box floating body
    int X120; //sphere floating bod
    int X131; //cylinder floating bod
    int X132; //cylinder floating bod
    int X133; //cylinder floating bod
    int X100; //delta x,y,z
    int X101; //ini Euler angles
    int X102; //ini linear velocity
    int X103; //ini angular velocity
    int X141;
    int X142;
    int X143;
    int X153; //symmetric wedge
    int X180; //read .stl file for floating body geometry
    int X181; //double scale .stl geometry
    int X182; //translation on/off
    int X183;
    int X210; //give fixed linear velocity
    int X211; //give fixed angular velocity
    int X310;
    int X311; //number of simple taut mooring lines
    int X312; //number of springs
    int X313; //initial rotation of mooring end points with 6DOF body
    int X314; //breaking mooring lines due to tension
    int X315; //breaking mooring lines due to time
    int X320;
    int X321; //number of nets
    int mooring_count;
    int net_count;
    double X21_d; //presribe homogeneous density floating body
    double X22_m; //prescribe mass floating body
    double X23_x;
    double X23_y;
    double X23_z;
    double X24_Ix;
    double X24_Iy;
    double X24_Iz;
    double X25_Cp; //damping rotation
    double X25_Cq; //damping rotation
    double X25_Cr; //damping rotation
    double X26_Cu; //damping translationa
    double X26_Cv; //damping translationa
    double X26_Cw; //damping translationa
    double X41; //eps for continuous forcing heavisi
    double X42; //distance for pressure force evaluation
    double X43; //distance for shear stress evaluation
    double X44; //viscosity in body
    double X100_x;
    double X100_y;
    double X100_z;
    double X101_phi;
    double X101_theta;
    double X101_psi;
    double X102_u;
    double X102_v;
    double X102_w;
    double X103_p;
    double X103_q;
    double X103_r;
    double *X110_xs;
    double *X110_xe;
    double *X110_ys;
    double *X110_ye;
    double *X110_zs;
    double *X110_ze;
    double X120_rad;
    double X120_xc;
    double X120_yc;
    double X120_zc;
    double X131_rad;
    double X131_h;
    double X131_xc;
    double X131_yc;
    double X131_zc;
    double X132_rad;
    double X132_h;
    double X132_xc;
    double X132_yc;
    double X132_zc;
    double X133_rad;
    double X133_h;
    double X133_xc;
    double X133_yc;
    double X133_zc;
    double X153_xs;
    double X153_xe;
    double X153_ys;
    double X153_ye;
    double X153_zs;
    double X153_ze;
    int X163; //wedge
    double *X163_x1;
    double *X163_y1;
    double *X163_z1;
    double *X163_x2;
    double *X163_y2;
    double *X163_z2;
    double *X163_x3;
    double *X163_y3;
    double *X163_z3;
    double *X163_x4;
    double *X163_y4;
    double *X163_z4;
    double *X163_x5;
    double *X163_y5;
    double *X163_z5;
    double *X163_x6;
    double *X163_y6;
    double *X163_z6;
    int X164; //hexahedro
    double *X164_x1;
    double *X164_y1;
    double *X164_z1;
    double *X164_x2;
    double *X164_y2;
    double *X164_z2;
    double *X164_x3;
    double *X164_y3;
    double *X164_z3;
    double *X164_x4;
    double *X164_y4;
    double *X164_z4;
    double *X164_x5;
    double *X164_y5;
    double *X164_z5;
    double *X164_x6;
    double *X164_y6;
    double *X164_z6;
    double *X164_x7;
    double *X164_y7;
    double *X164_z7;
    double *X164_x8;
    double *X164_y8;
    double *X164_z8;
    double X181_x; //scaling of stl geometry
    double X181_y; //scaling of stl geometry
    double X181_z; //scaling of stl geometry
    double X182_x; //translation of stl geometry
    double X182_y; //translation of stl geometry
    double X182_z; //translation of stl geometry
    double X183_x;
    double X183_y;
    double X183_z;
    double X183_phi;
    double X183_theta;
    double X183_psi;
    int X185; //stl refineme
    int X188; //ray cast algorithm
    double X186; //refinement factor
    int X205; //ype of ramp up function
    int X206; //ramp up velocity
    int X207; //ramp up draf
    double X206_ts; //ramp start
    double X206_te;
    double X207_ts; //ramp start
    double X207_te;
    double X210_u; //fixed u v
    double X210_v; //fixed v v
    double X210_w; //fixed w v
    double X211_p;
    double X211_q;
    double X211_r;
    int X240; //read 6DOF motion file
    double X241; //delta t for motion fi
    double X242_x; //delta x for motion fi
    double X242_y; //delta x for motion fi
    double X242_z; //delta x for motion fi
    double X243; //delta CoG for motion fi
    double *X311_xs;
    double *X311_xe;
    double *X311_ys;
    double *X311_ye;
    double *X311_zs;
    double *X311_ze;
    double *X311_w;
    double *X311_rho_c;
    double *X311_EA;
    double *X311_d;
    double *X311_l;
    double *X311_H;
    double *X311_P;
    double *X311_facT;
    double *X312_k;
    double *X312_T0;
    double *X314_T;
    double *X315_t;
    int *X320_type;
    double *X321_Sn;
    double *X321_d;
    double *X321_lambda;
    double *X321_dk;
    double *X321_rho;
    double *X321_nd;
    double *X321_nl;
    double *X322_D;
    double *X322_L;
    double *X322_x0;
    double *X322_y0;
    double *X322_z0;
    double *X322_phi;
    double *X322_theta;
    double *X322_psi;
    int X324;
    double X323_m; //dynamic net sinker properties
    double X323_d; //dynamic net sinker properties
    double X323_l; //dynamic net sinker properties
    double *X324_x;
    double *X324_y;
    double *X324_z;
    double X325_dt; //dynamic net time step
    double X325_relX; //dynamic net relaxation factors
    double X325_relY; //dynamic net relaxation factors
    double X325_relZ; //dynamic net relaxation factors
    int X400; //sflow external pressure term
    double X401_p0; //sflow external pressure term p0
    double X401_cl; //sflow external pressure term c
    double X401_cb; //sflow external pressure term c
    double X401_a; //sflow external pressure term a

    // FSI
    int Z10; //turn FSI o
    int Z11;
    int FSI_count;
    double *Z11_x;
    double *Z11_y;
    double *Z11_z;
    double *Z11_l;
    double *Z11_w;
    double *Z11_t;
    double *Z11_rho;
    double *Z11_e;
    double *Z11_ix;
    double *Z11_iy;
    double *Z11_iz;
    double *Z11_nu;
    double *Z11_n;
    double Z12_ckx; //fsi beam structural damping coefficients
    double Z12_cky; //fsi beam structural damping coefficients
    double Z12_ckz; //fsi beam structural damping coefficients
    double Z12_cdx; //fsi beam structural damping coefficients
    double Z12_cdy; //fsi beam structural damping coefficients
    double Z12_cdz; //fsi beam structural damping coefficients
	
	// Grid
    int Y40;
    int Y50;
    int Y60; //require
    int Y71; //turn on/off solid gcparax
    int Y72; //turn on/off solid gcparax
    int Y73; //turn on/off solid gcparax
    int Y74; //turn on/off solid gcparax

    // Test options
    int Y1; //turn on/off experimental screen force model
    int Y2; //turn external moments on/off
    int Y3;
    int Y4;
    int Y5;

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
