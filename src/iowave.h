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

#ifndef IOWAVE_RELAX_H_
#define IOWAVE_RELAX_H_

#include"ioflow.h"
#include"wave_interface.h"
#include"increment.h"
#include"flowfile_in.h"
#include"field4.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"

class fdm_fnpf;
class patchBC_interface;
class linear_regression_cont;

using namespace std;

class iowave final : public ioflow, public wave_interface, public increment, public flowfile_in
{

public:
    iowave(lexer*, ghostcell*,patchBC_interface*);
    virtual ~iowave();
    void gcio_update(lexer*,fdm*,ghostcell*) override final;
    void gcio_update_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void inflow_walldist(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*) override final {};
    void fsfinflow(lexer*,fdm*,ghostcell*) override final;
    void discharge(lexer*,fdm*,ghostcell*) override final;
    void inflow(lexer*,fdm*,ghostcell*,field&,field&,field&) override final;
    void inflow_plain(lexer*,fdm*,ghostcell*,field&,field&,field&);
    void rkinflow(lexer*,fdm*,ghostcell*,field&,field&,field&) override final;
    void fsfrkin(lexer*,fdm*,ghostcell*,field&) override final;
    void fsfrkout(lexer*,fdm*,ghostcell*,field&) override final;
    void isource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void jsource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void ksource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void pressure_io(lexer*,fdm*,ghostcell*) override final;
    void turbulence_io(lexer*,fdm*,ghostcell*) override final;
    void veltimesave(lexer*,fdm*,ghostcell*,vrans*) override final;

    void flowfile(lexer*,fdm*,ghostcell*,turbulence*) override final;

    void u_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void v_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void w_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void p_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void phi_relax(lexer*,ghostcell*,field&) override final;
    void vof_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void turb_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void U_relax(lexer*,ghostcell*,double*,double*) override final;
    void V_relax(lexer*,ghostcell*,double*,double*) override final;
    void W_relax(lexer*,ghostcell*,double*,double*) override final;
    void P_relax(lexer*,ghostcell*,double*) override final;
    void WL_relax(lexer*,ghostcell*,slice&,slice&) override final;
    void fivec_relax(lexer*, ghostcell*, double*) override final;
    void fifsf_relax(lexer*, ghostcell*, slice&) override final;
    void test_relax(lexer*, ghostcell*, slice&)override final;
    void eta_relax(lexer*,ghostcell*,slice&) override final;
    void um_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void vm_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void wm_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void ws_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void pm_relax(lexer*,ghostcell*,slice&) override final;

    // 2D
    void wavegen_2D_precalc(lexer*,fdm2D*,ghostcell*) override final;
    void wavegen_2D_precalc_ini(lexer*,ghostcell*) override final;

    void discharge2D(lexer*,fdm2D*,ghostcell*) override final;
    void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&) override final {};
    void Qin2D(lexer*,fdm2D*,ghostcell*) override final {};
    void Qout2D(lexer*,fdm2D*,ghostcell*) override final {};
    void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override final;
    void rkinflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override final;
    void isource2D(lexer*,fdm2D*,ghostcell*) override final;
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final;
    void full_initialize2D(lexer*,fdm2D*,ghostcell*) override final;

    double wave_fsf(lexer*,ghostcell*,double) override final;
    double wave_xvel(lexer*,ghostcell*,double,double,double) override final;
    double wave_yvel(lexer*,ghostcell*,double,double,double) override final;
    double wave_zvel(lexer*,ghostcell*,double,double,double) override final;

    int iozonecheck(lexer*,fdm*) override final;

    void ini(lexer*,fdm*,ghostcell*) override final;
    void ini_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void ini_fnpf(lexer*,fdm_fnpf*,ghostcell*) override final;
    void ini2D(lexer*,fdm2D*,ghostcell*) override final;
    void ini_ptf(lexer*,fdm*,ghostcell*) override final;

    void vrans_sed_update(lexer*,fdm*,ghostcell*,vrans*) override final;

    void waterlevel_update(lexer*,fdm*,ghostcell*) override final;

    // precalc
    void wavegen_precalc(lexer*,ghostcell*) override final;
    void wavegen_precalc_ini(lexer*,ghostcell*) override final;

    // FNPF
    void wavegen_precalc_fnpf(lexer*,fdm_fnpf*,ghostcell*) override final;
    void inflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,double*,double*,slice&,slice&) override final;
    void rkinflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;

    // NHFLOW
    void wavegen_precalc_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void wavegen_precalc_ini_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override final;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override final;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*) override final;
    void isource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void jsource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void ksource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&) override final;
    void turb_relax_nhflow(lexer*,fdm_nhf*,ghostcell*,double*) override final;

private:
    void Qin(lexer*,fdm*,ghostcell*);
    void Qout(lexer*,fdm*,ghostcell*);

    void gcio_update_impl(lexer*);

    void hydrograph_in_read(lexer*,fdm*,ghostcell*);
    void hydrograph_out_read(lexer*,fdm*,ghostcell*);
    double hydrograph_ipol(lexer*,ghostcell*,double**,int);

    void wavegen_precalc_space(lexer*,ghostcell*);
    void wavegen_precalc_space_dirichlet(lexer*,ghostcell*);
    void wavegen_precalc_time(lexer*,ghostcell*);
    void wavegen_precalc_decomp_relax(lexer*,ghostcell*);

    void wavegen2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&);
    void active_beach2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&);
    void inflow2D_plain(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&);

    void full_initialize(lexer*,fdm*,ghostcell*);
    void full_initialize_fnpf(lexer*,fdm_fnpf*,ghostcell*);
    void full_initialize_ptf(lexer*,fdm*,ghostcell*);
    void active_beach(lexer*,fdm*,ghostcell*,field&,field&,field&);
    void active_wavegen(lexer*,fdm*,ghostcell*,field&,field&,field&);
    void dirichlet_wavegen(lexer*,fdm*,ghostcell*,field&,field&,field&);

    void velini(lexer*,fdm*,ghostcell*);
    void pressure_outlet(lexer*,fdm*,ghostcell*);
    void pressure_inlet(lexer*,fdm*,ghostcell*);
    double local_fsf(lexer*,fdm*,ghostcell*);

    void gen_ini(lexer*,fdm*,ghostcell*);

    void wavegen_precalc_relax(lexer*,ghostcell*);
    void wavegen_precalc_relax_ini(lexer*,ghostcell*);
    void wavegen_precalc_dirichlet(lexer*,ghostcell*);
    void wavegen_precalc_dirichlet_ini(lexer*,ghostcell*);

    void wavegen_precalc_relax_func(lexer*,ghostcell*);
    void wavegen_precalc_relax_func_fnpf(lexer*,ghostcell*);
    void wavegen_precalc_relax_func_nhflow(lexer*,ghostcell*);

    // FNPF
    void fnpf_precalc_relax(lexer*,ghostcell*);
    void fnpf_precalc_relax_ini(lexer*,ghostcell*);
    void fnpf_precalc_parallel_relax(lexer*,ghostcell*);
    void fnpf_precalc_parallel_relax_ini(lexer*,ghostcell*);
    void fnpf_precalc_dirichlet(lexer*,ghostcell*);
    void fnpf_precalc_dirichlet_ini(lexer*,ghostcell*);
    void dirichlet_wavegen_fnpf(lexer*,fdm_fnpf*,ghostcell*,double*,double*,slice&,slice&);
    void active_beach_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, double*, slice&, slice&);

    void wavegen_precalc_decomp_space_fnpf(lexer*,ghostcell*);
    void wavegen_precalc_decomp_space_dirichlet_fnpf(lexer*,ghostcell*);
    void wavegen_precalc_decomp_time_fnpf(lexer*,ghostcell*);
    void wavegen_precalc_decomp_relax_fnpf(lexer*,ghostcell*);
    void wavegen_precalc_decomp_dirichlet_fnpf(lexer*,ghostcell*);

    void nhflow_precalc_relax(lexer*,fdm_nhf*,ghostcell*);
    void nhflow_precalc_relax_ini(lexer*,fdm_nhf*,ghostcell*);
    void nhflow_precalc_dirichlet(lexer*,fdm_nhf*,ghostcell*);
    void nhflow_precalc_dirichlet_ini(lexer*,fdm_nhf*,ghostcell*);

    void nhflow_dirichlet_wavegen(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&);
    void nhflow_active_wavegen(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&);
    void nhflow_active_beach(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*);
    void nhflow_inflow_plain(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*);
    void full_initialize_nhflow(lexer*,fdm_nhf*,ghostcell*);

    void nhflow_wavegen_precalc_decomp_space(lexer*,ghostcell*);
    void nhflow_wavegen_precalc_decomp_space_dirichlet(lexer*,ghostcell*);
    void nhflow_wavegen_precalc_decomp_time(lexer*,ghostcell*);
    void nhflow_wavegen_precalc_decomp_relax(lexer*,fdm_nhf*,ghostcell*);
    void nhflow_wavegen_precalc_decomp_dirichlet(lexer*,ghostcell*);

    void Qin_nhf(lexer*,fdm_nhf*,ghostcell*);
    void Qout_nhf(lexer*,fdm_nhf*,ghostcell*);

    void timeseries(lexer*,ghostcell*);

    slice4 eta;
    // Depth-dependent: wave_u/v/w and the level-set distance all take z,
    // so these are 3D fields staggered like the velocity they relax.
    field1 uval_amrex;
    field2 vval_amrex;
    field3 wval_amrex;
    field4 lsval_amrex;

    slice1 relax1_wg, relax1_nb;
    slice2 relax2_wg, relax2_nb;
    slice4 relax4_wg, relax4_nb;
    sliceint4 wgflag;

    double rb1(lexer*,double);
    double rb3(lexer*,double);

    double rb1_ext(lexer*,int);
    double rb3_ext(lexer*,int);
    
    int rb1_flag(lexer*,int);

    double ramp(lexer*);
    double ramp_corr(lexer*);

    double xgen(lexer*);
    double xgen1(lexer*);
    double xgen2(lexer*);
    double ygen(lexer*);
    double ygen1(lexer*);
    double ygen2(lexer*);

    double distgen(lexer*);
    double distbeach(lexer*);

    void distbeach_ini(lexer*);
    void distgen_ini(lexer*);
    int intriangle(lexer*,double,double,double,double,double,double,double,double);

    //PLIC
    double V0Calc_PLIC(lexer*, fdm*, double, double, double, double);
    slice4 vofheight;
    slice4 genheight;
    field4 vofgen;

    int n,count;
    int wtype;
    double inflow_bed,uvel,vvel,wvel;
    double uhvel,vhvel,whvel;
    double area,Ai,Ao,Ui,fac;
    double dist1,dist2,dist2_fac;
    double x,y,z;
    double x1,y1,x2,y2,z3;
    double xg,yg,zg,dg,db;
    double xc,yc,zc;
    int gcval_press;
    const double epsi,psi;
    double alpha,*beta,gamma;
    double H,G,phival;
    double kinval,epsval;
    double tan_alpha,*tan_beta;
    double wh;
    int beach_relax;
    double starttime;

    int gcgen1_count,gcgen2_count,gcgen3_count,gcgen4_count;
    int **gcgen1,**gcgen2,**gcgen3,**gcgen4;

    // relax points
    double **G1,**G2,**G3,**G4,**Gs,**Ge;
    double **B1,**B2,**B3,**B4,**Bs,**Be;

    // relax pre-calc
    int wave_comp;
    int upt_count,vpt_count,wpt_count,ppt_count,ept_count;
    std::vector<double> uval, vval, wval, lsval, Fival, Fifsfval, Fifsfval0, Fifsfval1, Uinval;
    std::vector<double> UHval, VHval, WHval;
    std::vector<double> vofval;

    // decomp space
    std::vector<std::vector<double>> etaval_S_sin,Fifsfval_S_sin;
    std::vector<std::vector<float>> uval_S_sin, vval_S_sin, wval_S_sin, Fival_S_sin;
    std::vector<std::vector<double>> etaval_S_cos,Fifsfval_S_cos;
    std::vector<std::vector<float>> uval_S_cos, vval_S_cos, wval_S_cos, Fival_S_cos;
    
    // decomp time
    std::vector<double> uval_T_sin, vval_T_sin, wval_T_sin, etaval_T_sin, Fival_T_sin, Fifsfval_T_sin;
    std::vector<double> uval_T_cos, vval_T_cos, wval_T_cos, etaval_T_cos, Fival_T_cos, Fifsfval_T_cos;

    double zloc1,zloc2,zloc3,zloc4,zcoor;

    slice4 wsfmax;
    double time_n,time_0,time_1;

    bool u_switch,v_switch,w_switch,p_switch,h_switch,f_switch;

    double **hydro_in,**hydro_out;
    int hydro_in_count,hydro_out_count;

    double netQ,netQ_n,netV;
    double netV_corr,netV_corr_n;
    double b0,b1;
    double timemeter;

    linear_regression_cont *linreg;
    patchBC_interface *pBC;
    
    double cosh_func(double);
    
};

#endif


