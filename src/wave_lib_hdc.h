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

#ifndef WAVE_LIB_HDC_H_
#define WAVE_LIB_HDC_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"

#include<fstream>
#include<vector>

using namespace std;

class wave_lib_hdc final : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    wave_lib_hdc(lexer*, ghostcell*);
    virtual ~wave_lib_hdc();

    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;

    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;

private:
    // functions
    void wave_prestep_cfd(lexer*,ghostcell*);
    void wave_prestep_fnpf(lexer*,ghostcell*);

    void read_header(lexer*,ghostcell*);
    void read_result_cfd(lexer*,ghostcell*,std::vector<std::vector<double>>&,std::vector<std::vector<std::vector<double>>>&,std::vector<std::vector<std::vector<double>>>&,std::vector<std::vector<std::vector<double>>>&,int);
    void read_result_fnpf(lexer*,ghostcell*,std::vector<std::vector<double>>&,std::vector<std::vector<double>>&,int);

    void fill_conti_fnpf(lexer*,ghostcell*);
    void fill_conti_cfd(lexer*,ghostcell*);

    void filename_header(lexer*,ghostcell*);
    void filename_single(lexer*,ghostcell*,int);
    void filename_continuous(lexer*,ghostcell*);

    void allocate_fnpf();
    void allocate_cfd();

    void time_interpol_fnpf(lexer*);
    void time_interpol_cfd(lexer*);

    // interpolation
    double ccpol3D(lexer*,std::vector<std::vector<std::vector<double>>>&,double,double,double);
    double ccpol2D(lexer*,std::vector<std::vector<double>>&,double,double);
    double ccpol2DM(lexer*,std::vector<std::vector<std::vector<double>>>&,double,double);
    double space_interpol(lexer*,std::vector<std::vector<std::vector<double>>>&,double,double,double);
    double plane_interpol(lexer*,std::vector<std::vector<double>>&,double,double);

    int pos_i(lexer*,double);
    int pos_j(lexer*,double);
    int pos_k(lexer*,double,int,int);

    int ihalf(int,int);

    // arrays
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Zsig;
    std::vector<std::vector<std::vector<double>>> Z;
    std::vector<std::vector<double>> B;
    std::vector<double> simtime;

    std::vector<std::vector<std::vector<double>>> U1, U2, U;
    std::vector<std::vector<std::vector<double>>> V1, V2, V;
    std::vector<std::vector<std::vector<double>>> W1, W2, W;
    std::vector<std::vector<double>> E1, E2, E;
    std::vector<std::vector<double>> F1, F2, F;

    // variables
    double t_start,t_end;
    int endseries;
    double val;
    int q1,q2,q1n,q2n;
    double t1,t2,tn,deltaT;
    int Nx,Ny,Nz;
    int num;

    int iin;
    float ffn;
    double ddn;
    char name[300];
    ifstream result;

    int file_version,file_type,file_conti;
    int ip1,jp1,kp1;
    int ii,jj,kk;
    int iii,jjj,kkk;
    double xp,yp,zp;
    double wa,wb,wc;
    double v1,v2,v3,v4,v5,v6,v7,v8;
    double x1,x2,x3,x4,y1,y2;

    double Xstart,Xend,Ystart,Yend;

    int startup;
    int numiter,diter,jdir;
    int file_iter;
};

#endif
