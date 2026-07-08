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

#include "driver.h"
#include "lexer.h"
#include "ghostcell.h"
#include "freesurface_header.h"
#include "turbulence_header.h"
#include "momentum_header.h"
#include "pressure_header.h"
#include "fdm_header.h"
#include "sediment_header.h"
#include "convection_header.h"
#include "solver_header.h"
#include "field_header.h"
#include "heat_header.h"
#include "concentration_header.h"
#include "benchmark_header.h"
#include "6DOF_header.h"
#include "FSI_header.h"
#include "vrans_header.h"
#include "waves_header.h"

void driver::logic_cfd()
{
    makegrid_cds();
    pini = new initialize(p);

    if(p->mpirank==0)
    cout<<"starting ini"<<endl;
    pini->start(p,a,pgc);

    if(p->mpirank==0)
    cout<<"creating objects"<<endl;

    // time stepping
    switch (p->N48)
    {
    case 0:
        ptstep = new fixtimestep(p);
        break;

    case 1:
        if(p->D20==1)
        ptstep = new etimestep(p);
        else
        ptstep = new ietimestep(p);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid time stepping method (N 48) specified in input file: " << p->N48 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Multiphase
    if(p->F300>0)
    pmp = new multiphase_f(p,a,pgc);
    else
    pmp = new multiphase_v();

    // discretization scheme
    // Convection
    switch (p->D10)
    {
    case 0:
        pconvec = new convection_void(p);
        break;
    case 1:
        pconvec = new fou(p);
        break;
    case 2:
        pconvec = new cds2(p);
        break;
    case 3:
        pconvec = new quick(p);
        break;
    case 4:
        pconvec = new weno_flux_nug(p);
        break;
    case 5:
        pconvec = new weno_hj_nug(p);
        break;
    case 6:
        pconvec = new cds4(p);
        break;
    case 7:
        pconvec = new weno3_flux(p);
        break;
    case 8:
        pconvec = new weno3_hj(p);
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
    case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
        pconvec = new hires(p,p->D10);
        break;
    case 60:
        pconvec = new hcds6(p);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid convection scheme (D 10) specified in input file: " << p->D10 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Convection Turbulence
    switch (p->T12)
    {
    case 0:
        pturbdisc = new convection_void(p);
        break;
    case 1:
        pturbdisc = new ifou(p);
        break;
    case 5:
        pturbdisc = new iweno_hj_df_nug(p);
        break;
    case 55:
        pturbdisc = new iweno_hj(p);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid turbulence convection scheme (T 12) specified in input file: " << p->T12 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Convection FSF
    if(p->F85==0)
    {
        switch (p->F35)
        {
        case 0:
            pfsfdisc = new convection_void(p);
            break;
        case 1:
            pfsfdisc = new fou(p);
            break;
        case 2:
            pfsfdisc = new cds2_alt(p);
            break;
        case 3:
            pfsfdisc = new quick(p);
            break;
        case 4:
            pfsfdisc = new weno_flux_nug(p);
            break;
        case 5:
            pfsfdisc = new weno_hj_nug(p);
            break;
        case 6:
            pfsfdisc = new cds4(p);
            break;
        case 7:
            pfsfdisc = new weno3_flux(p);
            break;
        case 8:
            pfsfdisc = new weno3_hj(p);
            break;
        case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
        case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
            pfsfdisc = new hires(p,p->F35);
            break;
        case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
            pfsfdisc = new hires(p,p->F35);
            break;
        default:
            if(p->mpirank==0)
            std::cerr << "Invalid convection scheme for free surface (F 35) specified in input file: " << p->F35 << std::endl;
            pgc->final(EXIT_FAILURE);
        }
    }
    // Convection VOF
    else
    {
        switch (p->F85)
        {
        case 1:
            pfsfdisc = new fou(p);
            break;
        case 2:
            pfsfdisc = new cds2_alt(p);
            break;
        case 3:
            pfsfdisc = new quick(p);
            break;
        case 4:
            pfsfdisc = new weno_flux_nug(p);
            break;
        case 5:
            pfsfdisc = new weno_hj_nug(p);
            break;
        case 6:
            pfsfdisc = new cds4(p);
            break;
        case 7:
            pfsfdisc = new weno3_flux(p);
            break;
        case 8:
            pfsfdisc = new weno3_hj(p);
            break;
        case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
        case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
            pfsfdisc = new hires(p,p->F85);
            break;
        case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
            pfsfdisc = new hires(p,p->F85);
            break;
        case 51:
            pfsfdisc = new hric(p);
            break;
        case 52:
            pfsfdisc = new hric_mod(p);
            break;
        case 53:
            pfsfdisc = new cicsam(p);
            break;
        default:
            if(p->mpirank==0)
            std::cerr << "Invalid convection scheme for VOF (F 85) specified in input file: " << p->F85 << std::endl;
            pgc->final(EXIT_FAILURE);
        }
    }

    // Convection Multiphase LSM
    switch (p->F305)
    {
    case 0:
        pmpconvec = new convection_void(p);
        break;
    case 1:
        pmpconvec = new fou(p);
        break;
    case 2:
        pmpconvec = new cds2(p);
        break;
    case 3:
        pmpconvec = new quick(p);
        break;
    case 4:
        pmpconvec = new weno_flux_nug(p);
        break;
    case 5:
        pmpconvec = new weno_hj_nug(p);
        break;
    case 6:
        pmpconvec = new cds4(p);
        break;
    case 7:
        pmpconvec = new weno3_flux(p);
        break;
    case 8:
        pmpconvec = new weno3_hj(p);
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
    case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
        pmpconvec = new hires(p,p->F305);
        break;
    case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
        pmpconvec = new hires(p,p->F305);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid convection scheme for multiphase level set (F 305) specified in input file: " << p->F305 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Convection Concentration

    if(p->S12>=1)
        pconcdisc = new iweno_hj(p);
    else
    {
        switch (p->C15)
        {
        case 0:
            pconcdisc = new convection_void(p);
            break;
        case 1:
            pconcdisc = new fou(p);
            break;
        case 2:
            pconcdisc = new cds2_alt(p);
            break;
        case 3:
            pconcdisc = new quick(p);
            break;
        case 4:
            pconcdisc = new weno_flux_nug(p);
            break;
        case 5:
            pconcdisc = new weno_hj_nug(p);
            break;
        case 6:
            pconcdisc = new cds4(p);
            break;
        case 7:
            pconcdisc = new weno3_flux(p);
            break;
        case 8:
            pconcdisc = new weno3_hj(p);
            break;
        case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
        case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
            pconcdisc = new hires(p,p->C15);
            break;
        case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
            pconcdisc = new hires(p,p->C15);
            break;
        default:
            if(p->mpirank==0)
            std::cerr << "Invalid convection scheme (C 15) specified in input file: " << p->C15 << std::endl;
            pgc->final(EXIT_FAILURE);
        }
    }

    // turbulence model
    switch (p->T10)
    {
    case 0:
        pturb = new kepsilon_void();
        break;
    case 1: case 21:
        pturb = new kepsilon_IM1(p,a,pgc);
        break;
    case 2: case 22:
        if(p->F80==4)
        pturb = new komega_IM1_PLIC(p,a,pgc);
        else
        pturb = new komega_IM1(p,a,pgc);
        break;
    case 12:
        pturb = new EARSM_kw_IM1(p,a,pgc);
        break;
    case 31:
        pturb = new LES_smagorinsky(p,a);
        break;
    case 33:
        pturb = new LES_WALE(p,a);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid turbulence model (T 10) specified in input file: " << p->T10 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Heat
    switch (p->H10)
    {
    case 0:
        pheat = new heat_void(p,a,pgc);
        break;
    case 1:
        pheat = new heat_AB(p,a,pgc,pheat);
        break;
    case 2:
        pheat = new heat_RK2(p,a,pgc,pheat);
        break;
    case 3:
        pheat = new heat_RK3(p,a,pgc,pheat);
        break;
    case 4:
        pheat = new heat_RK3CN(p,a,pgc,pheat);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid heat equation time discretization scheme (H 10) specified in input file: " << p->H10 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Convection Heat
    switch (p->H15)
    {
    case 0:
        pheatdisc = new convection_void(p);
        break;
    case 1:
        pheatdisc = new fou(p);
        break;
    case 2:
        pheatdisc = new cds2(p);
        break;
    case 3:
        pheatdisc = new quick(p);
        break;
    case 4:
        pheatdisc = new weno_flux_nug(p);
        break;
    case 5:
        pheatdisc = new weno_hj_nug(p);
        break;
    case 6:
        pheatdisc = new cds4(p);
        break;
    case 7:
        pheatdisc = new weno3_flux(p);
        break;
    case 8:
        pheatdisc = new weno3_hj(p);
        break;
    case 9:
        pheatdisc = new weno_flux(p);
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
    case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
        pheatdisc = new hires(p,p->H15);
        break;
    case 60:
        pheatdisc = new hcds6(p);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid convection scheme for heat equation (H 15) specified in input file: " << p->H15 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Concentration
    switch (p->C10)
    {
    case 0:
        pconc = new concentration_void(p,a,pgc);
        break;
    case 1:
        pconc = new concentration_AB(p,a,pgc);
        break;
    case 2:
        pconc = new concentration_RK2(p,a,pgc);
        break;
    case 3:
        pconc = new concentration_RK3(p,a,pgc);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid concentration equation time discretization scheme (C 10) specified in input file: " << p->C10 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

// Diffusion
    // momentum and scalars
    switch (p->D20)
    {
    case 0:
        pdiff = new diff_void();
        break;
    case 1:
        pdiff = new ediff2(p);
        break;
    case 2:
        if(p->F80==4)
        {
            if(p->j_dir)
            pdiff = new idiff2_PLIC(p);
            else
            pdiff = new idiff2_PLIC_2D(p);
        }
        else
        {
            if(p->j_dir)
            pdiff = new idiff2_FS(p);
            else
            pdiff = new idiff2_FS_2D(p);
        }
        break;
    case 3:
        if(p->j_dir)
        pdiff = new idiff2_CN(p);
        else
        {
            if(p->mpirank==0)
            std::cerr << "Diffusion scheme (D 20) specified in input file: " << p->D20 << " is not implemented for 2D. Please set D 20 to 0, 1, or 2 for 2D simulations." << std::endl;
            pgc->final(EXIT_FAILURE);
        }
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid diffusion scheme (D 20) specified in input file: " << p->D20 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // turbulence
    if(p->T10>0 && p->D20>0)
    pturbdiff = new idiff2(p);
    else
    pturbdiff = new diff_void;

    // concentration
    if(p->D20==0 || p->C10==0)
    pconcdiff = new diff_void;
    else if(p->D20==1 && p->C10<=10 && p->C10>0)
    pconcdiff = new ediff2(p);
    else if(p->D20>=2 && p->C10<=10 && p->C10>0)
    pconcdiff = new idiff2_FS(p);

    // Free Surface
    // LSM
    if((p->F30==0 || (p->N40==2 || p->N40==3 || p->N40==4 || p->N40==22 || p->N40==23 || p->N40==24 || p->N40==33)) && p->F80==0)
    pfsf = new levelset_void(p,a,pgc,pheat,pconc);
    else if(p->F30==2 && (p->N40==44 || p->N40==12 || p->N40==13 || p->N40==14))
    pfsf = new levelset_RK2(p,a,pgc,pheat,pconc);
    else if(p->F30==3 && (p->N40==44 || p->N40==12 || p->N40==13 || p->N40==14))
    pfsf = new levelset_RK3(p,a,pgc,pheat,pconc);
    // VOF
    if(p->F80>0)
    {
        if(p->N40==2 || p->N40==3 || p->N40==22 || p->N40==23 || p->N40==33)
        {
            if(p->mpirank==0)
            std::cerr << "Invalid free surface advection method (N 40) specified in input file: " << p->N40 << " is not compatible with VOF convection scheme. Please set N 40 to 12, 13, or 14 for VOF simulations." << std::endl;
            pgc->final(EXIT_FAILURE);
        }

        if(p->F80==1)
        pfsf = new VOF_AB(p,a,pgc,pheat);
        else if(p->F80==3)
        pfsf = new VOF_RK3(p,a,pgc,pheat);
        else if(p->F80==4 && (p->N40==12 || p->N40==13))
        pfsf = new VOF_PLIC(p,a,pgc,pheat);
        else
        {
            if(p->mpirank==0)
            std::cerr << "Invalid VOF convection scheme (F 80) specified in input file: " << p->F80 << std::endl;
            pgc->final(EXIT_FAILURE);
        }
    }

    // Reini
    switch (p->F40)
    {
    case 0:
        preini = new reini_void(p);
        break;
    case 3: case 23:
        preini = new reini_RK3(p,1);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid reinitialization method (F 40) specified in input file: " << p->F40 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // pressure scheme
    switch (p->D30)
    {
    case 0:
        ppress = new pressure_void(p);
        break;
    case 1: case 2: case 3:
        ppress = new pjm_corr(p,a,pgc,pheat,pconc);
        // poisson scheme for pressure
        ppois = new poisson_pcorr(p,pheat,pconc);
        break;
    case 10:
        ppress = new pjm_hydrostatic(p,a,pheat,pconc);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid pressure scheme (D 30) specified in input file: " << p->D30 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Solver
    if(p->j_dir)
    {
        if(p->N10==40 || p->N10==41)
        psolv = new hypre_ssamg(p,a,pgc);
        else
        psolv = new bicgstab_ijk(p,a,pgc);
    }
    else
    {
        if(p->N10==40 || p->N10==41)
        psolv = new hypre_ssamg(p,a,pgc);
        else
        psolv = new bicgstab_ijk_2D(p,a,pgc);
    }

    // Poison Solver
    switch (p->N10)
    {
    case 0:
        ppoissonsolv = new solver_void(p,a,pgc);
        break;
    case 1:
        if(p->j_dir)
        ppoissonsolv = new bicgstab_ijk(p,a,pgc);
        else
        ppoissonsolv = new bicgstab_ijk_2D(p,a,pgc);
        break;
    case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19:
        ppoissonsolv = new hypre_struct(p,pgc,p->N10,p->N11);
        break;
    case 20: case 21: case 22: case 23: case 24: case 25: case 26: case 27: case 28: case 29:
        ppoissonsolv = new hypre_aij(p,a,pgc);
        break;
    case 30: case 31: case 32: case 33: case 34: case 35: case 36: case 37: case 38: case 39:
        ppoissonsolv = new hypre_sstruct(p,a,pgc);
        break;
    case 40: case 41:
        ppoissonsolv = new hypre_ssamg(p,a,pgc);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid Poisson solver (N 10) specified in input file: " << p->N10 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // VRANS
    switch (p->B200)
    {
    case 0:
        pvrans = new vrans_v(p,pgc);
        break;
    case 1:
        pvrans = new vrans_f(p,pgc);
        break;
    case 2:
        pvrans = new vrans_veg(p,pgc);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid VRANS model (B 200) specified in input file: " << p->B200 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // IOFlow
    if(p->B60>=1)
    pflow = new ioflow_f(p,pgc,pBC);
    if(p->B90>=1)
    pflow = new iowave(p,pgc,pBC);
    if(p->B180==1 || p->B191==1 || p->B192==1)
    pflow = new ioflow_gravity(p,pgc,pBC);
    if(p->B60==0 && p->B90==0 && p->B180==0)
    pflow = new ioflow_v(p,pgc,pBC);
    if(pflow==nullptr)
    {
        if(p->mpirank==0)
        std::cerr << "Invalid flow boundary condition specified in input file. \nPlease set B 60, B 90, or B 180 to 1 for flow boundary conditions or all to 0." << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Potential Flow Solver
    switch (p->I11)
    {
    case 0:
        potflow = new potential_v;
        break;
    case 1:
        potflow = new potential_f(p);
        break;
    case 2:
        potflow = new potential_water(p);
        break;
    default:
        if(p->mpirank==0)    std::cerr << "Invalid potential flow solver specified in input file: " << p->I11 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Benchmark
    switch (p->F150)
    {
    case 0:
        pbench = new benchmark_void;
        break;
    case 1:
        pbench = new benchmark_vortex(p,a,pgc);
        break;
    case 2:
        pbench = new benchmark_disk(p,a);
        break;
    case 21:
        pbench = new benchmark_disk_yz(p,a);
        break;
    case 22:
        pbench = new benchmark_disk_xy(p,a);
        break;
    case 3:
        pbench = new benchmark_vortex3D(p,a);
        break;
    case 4:
        pbench = new benchmark_TaylorGreen(p,a);
        break;
    case 11:
        pbench = new benchmark_convection(p,a);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid benchmark case (F 150) specified in input file: " << p->F150 << std::endl;
        pgc->final(EXIT_FAILURE);
    }

    // Printer
    pprint = new printer_CFD(p,a,pgc);

    // Data
    if(p->P150>0)
    pdata = new expdata_f(p,a,pgc);
    else
    pdata = new expdata_void(p,a,pgc);

    // Sediment
    if(p->S10>0)
    {
        if(p->Q10==0)
        psed = new sediment_f(p,pgc,pturb,pBC);
        else if(p->Q10==1)
        psed = new sediment_part(p,a,pgc,pturb,pBC);
    }
    else
    psed = new sediment_void();

    if(p->S10>0 || p->G1==1 || p->toporead==1)
    {
        if(p->G40==0)
        preto = new reinitopo_void();
        else if(p->G40>0)
        preto = new reinitopo_RK3(p);
    }

    if(p->solidread==0 || p->G40==0)
    preso = new reinitopo_void();
    else if(p->solidread==1 && p->G40>0)
    preso = new reinisolid_RK3(p);

    // 6DOF
    if(p->X10==0)
    p6dof = new sixdof_void(p,pgc);
    else if(p->X10==1)
    p6dof = new sixdof_cfd(p,a,pgc);

    // FSI
    if(p->Z10==0)
    pfsi = new fsi_void(p,pgc);
    else if(p->Z10==1)
    pfsi = new fsi_strips(p,pgc);

    // Velocities
    switch (p->N40)
    {
    case 0:
        pmom = new momentum_void();
        break;
    case 2: case 22:
        pmom = new momentum_FC2(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        break;
    case 3: case 23:
        if(p->F80==4)
        pmom = new momentum_FC3_PLIC(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        else
        pmom = new momentum_FC3(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        break;
    case 4: case 24:
        if(p->F80==4)
        {
            if(p->mpirank==0)
            std::cerr << "Momentum equation time discretization scheme (N 40) specified in input file: " << p->N40 << " with free surface scheme (F 80) specified as VOF-PLIC is not implemented. Please set N 40 to 12, 13, or 33 for VOF-PLIC simulations." << std::endl;
            pgc->final(EXIT_FAILURE);
        }
        else
        pmom = new momentum_FCLS3(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        break;
    case 5:
        pmom = new momentum_RK3CN(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsi);
        break;
    case 12:
        pmom = new momentum_RK2(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsi);
        break;
    case 13:
        if(p->F80==4)
        pmom = new momentum_RK3_PLIC(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,pfsi);
        else
        pmom = new momentum_RK3(p,a,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsi);
        break;
    case 14:
        if(p->X10==0 && p->Z10==0)
        {
            pmom_sf = new momentum_RKLS3_sf(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
            pmom = new momentum_void();
        }
        else if(p->X10==1 || p->Z10>0)
        {
            pmom_df = new momentum_RKLS3_df(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);
            pmom = new momentum_void();
        }
        else
        {
            if(p->mpirank==0)
            std::cerr << "Invalid momentum equation time discretization scheme (N 40) specified in input file: " << p->N40 << " for FSI simulations. Please set N 40 to 12 or 13 for FSI simulations." << std::endl;
            pgc->final(EXIT_FAILURE);
        }
        break;
    case 33:
        if(p->F80==4)
        pmom = new momentum_FCC3_PLIC(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        else
        pmom = new momentum_FCC3(p,a,pgc,pconvec,pfsfdisc,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pheat,pconc,preini,pfsi);
        break;
    case 44:
        pmom = new momentum_RKLS3(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow,pfsi);
        break;
    default:
        if(p->mpirank==0)
        std::cerr << "Invalid momentum equation time discretization scheme (N 40) specified in input file: " << p->N40 << std::endl;
        pgc->final(EXIT_FAILURE);
    }
}

void driver::patchBC_logic()
{
    if((p->B440>0 || p->B441>0 || p->B442>0) && p->A10==2)
        pBC = new patchBC_2D(p,pgc);

    else if((p->B440>0 || p->B441>0 || p->B442>0) && p->A10==6)
        pBC = new patchBC(p,pgc);

    else
        pBC = new patchBC_void(p);
}
