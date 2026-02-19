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

#ifndef GHOSTCELL_H_
#define GHOSTCELL_H_

#include "increment.h"
#include <mpi.h>
#include <vector>

class convection;
class fdm;
class fdm2D;
class fdm_fnpf;
class fdm_nhf;
class field;
class fieldint;
class ioflow;
class lexer;
class reini;
class slice;
class sliceint;
class vec2D;

using namespace std;

class ghostcell final : public increment
{
public:
    ghostcell(int&,char**,lexer*);

    ghostcell(const ghostcell&) = delete;
    ghostcell& operator=(const ghostcell&) = delete;
    ghostcell(ghostcell&&) = delete;
    ghostcell& operator=(ghostcell&&) = delete;

    ~ghostcell() = default;

    void final(bool error=false);
    void gc_ini(lexer*);
    void gcx_ini(lexer*);
    void gcx_cart_topology(lexer*);
    void mpi_check(lexer*);

    void start1(lexer*,field&,int);
    void start2(lexer*,field&,int);
    void start3(lexer*,field&,int);
    void start4(lexer*,field&,int);
    void start4_sum(lexer*,field&,int);

    void start1V(lexer*,double*,int);
    void start2V(lexer*,double*,int);
    void start3V(lexer*,double*,int);
    void start4V(lexer*,double*,int);
    void start4V_par(lexer*,double*,int);
    void start5V(lexer*,double*,int);
    void start5Vfull(lexer*,double*,int);

    void start20V(lexer*,double*,int);
    void start24V(lexer*,double*,int);
    void start30V(lexer*,double*,int);
    void start49V(lexer*,double*,int);
    void start60V(lexer*,double*,int);

    void start7V(lexer*,double*,sliceint&,int);
    void start7S(lexer*,double*,int);
    void start7P(lexer*,double*,int);

    void startintV(lexer*,int*,int);
    void startintVF(lexer*,int*,int);

    void dgcpol1(lexer*,field&);
    void dgcpol2(lexer*,field&);
    void dgcpol3(lexer*,field&);
    void dgcpol4(lexer*,field&);

    // particle
    void parapls(lexer*,double**,double**,int*,int*);
    void gcpartnum(int[6],int[6]);
    void gcpartx(int[6],int[6],double*[6],double*[6]);

    //  Update
    void facenbx(lexer*,fieldint&,std::vector<int>&);
    void flagx(lexer*,int*);
    void flagx7(lexer*,int*);
    void rangex(lexer*,std::vector<int>&,int);
    void gcxupdate(lexer*);

    void rownum4_update(lexer*,fieldint&);
    void rownum7_update(lexer*,int*);

    void sizeM_update(lexer*);

    void fdm2D_update(fdm2D*);
    void fdm_fnpf_update(fdm_fnpf*);
    void fdm_nhf_update(fdm_nhf*);
    void fdm_update(fdm*);

    void gcb_velflagio(lexer*);

    // Forcing CFD
    void solid_forcing(lexer*,fdm*,double,field&,field&,field&,field&,field&,field&);
    void solid_forcing_ini(lexer*,fdm*);
    void solid_forcing_flag_update(lexer*,fdm*);
    void solid_forcing_lsm(lexer*,fdm*,field&);
    void solid_forcing_eta(lexer*,slice&);
    void solid_forcing_bed(lexer*,slice&);
    double Hsolidface(lexer*, fdm*, int,int,int);

    // 6DOF update gcdf
    void gcdf_update(lexer*,fdm*);
    void gcsldf_update(lexer*);

    // IBM
    void flagfield(lexer*);
    void flagbase(lexer*);

    // PARALLEL
    #if not USE_AMREX
    void gcparax(lexer*, field&, int);
    #endif
    void gcparaxint(lexer*, fieldint&, int);
    void gcparaxijk(lexer*, double*, int);
    void gcparaxijk_single(lexer*, double*, int);
    void gcparax7(lexer*, double*&, int);
    void gcparax7co(lexer*, double*, int);
    void gcparax7int(lexer*, int*&, int);
    void gcparax4a(lexer*, field&, int);
    void gcparax4_sum(lexer*, field&, int);
    void gcparacox4a_sum(lexer*, field&, int);
    void gcparaxV(lexer*, double*, int);
    void gcparaxintV(lexer*, int*, int);
    void gcparaxV1(lexer*, double*, int);
    void gcparacox(lexer*, field&);
    void gcparacoxV(lexer*, double*, int);
    void gcparacoxV1(lexer*, double*, int);
    void gcperiodicx(lexer*, field&, int);
    void gcsync();
    void verticalmax(lexer*,fdm*,slice&);
    double timer();

    //Collective Communication
    void gather_int(int *, int, int *, int);
    void gatherv_int(int*, int, int*, int*, int*);
    void allgather_int(int *, int, int *, int);
    void allgatherv_int(int *, int, int *, int*, int*);
    void gather_double(double *, int, double *, int);
    void gatherv_double(double *, int, double *, int*, int*);
    void bcast_int(int*, int);
    void bcast_double(double *, int, int=0);
    double globalsum(double);
    int globalisum(int);
    double globalmax(double);
    double globalmin(double);
    int globalimax(int);
    int globalimin(int);
    double timesync(double);
    void globalctrl(lexer*);
    //Utilities
    void walldistance(lexer*,fdm*,convection*,reini*,ioflow*,field&);

    MPI_Comm mpi_comm = MPI_COMM_NULL;

    // Slice
    // epol
    void gcsl_start1(lexer*,slice&,int);
    void gcsl_start2(lexer*,slice&,int);
    void gcsl_start4(lexer*,slice&,int);
    void gcsl_start4a(lexer*,slice&,int);

    void gcsl_start1int(lexer*,sliceint&,int);
    void gcsl_start2int(lexer*,sliceint&,int);
    void gcsl_start4int(lexer*,sliceint&,int);
    void gcsl_start4Vint(lexer*,int*, int);

    void gcsl_setbc1(lexer*);
    void gcsl_setbc2(lexer*);
    void gcsl_setbc4(lexer*);
    void gcsl_setbcio(lexer*);

    void dgcslini1(lexer*);
    void dgcslini2(lexer*);
    void dgcslini4(lexer*);

    // parallel
    void gcslparax_fh(lexer*, slice&);
    void gcslparax_int(lexer*, sliceint&);

    void gcslflagx(lexer*, int*);
    void gcslparaxijk(lexer*, double*);
    void gcslparaxijk_single(lexer*, double*);

    void fivec(lexer*,double*,sliceint&);
    void fivec2D(lexer*,double*,sliceint&);
    void fivec_vel(lexer*,double*,sliceint&);
    void fivec2D_vel(lexer*,double*,sliceint&);
    void gc_periodic(field&,int);

    //NHFLOW
    void gciobc_update(lexer*,fdm_nhf*);

    void dirichlet_ortho(field&,double,int);

private:
    void dgcslpol1(lexer*, slice&);
    void dgcslpol2(lexer*, slice&);
    void dgcslpol4(lexer*, slice&);

    void gcslparax(lexer*, slice&);
    void gcslparaxV_int(lexer*, int*);

    void gcslparacox(lexer*, slice&);
    void gcslparacoxV_int(lexer*, int*);
    void gcslparacox_int(lexer*, sliceint&);

    void Sendrecv_double(int,int,int,int,int,int);
    void Sendrecv_int(int,int,int,int,int,int);
    void Sendrecv_1D(const void*[6],int[6],void*[6],int[6],MPI_Datatype);
    void Sendrecv_2D(const void*[6],int[6],void*[6],int[6],MPI_Datatype);
    void Sendrecv_3D(const void*[6],int[6],void*[6],int[6],MPI_Datatype);

    void set_DF(lexer *p, fdm *a);

    void gcwait(lexer*);

    // boundary conditions
    void heatbc(field&,int);
    void neumann(field&,int);
    void noslip(field&,int);
    void outflow(field&,int);
    void potentialbc(field&,int);

    // Slice BCs
    void gcsl_neumann(slice&,int);
    void gcsl_neumann_hx(slice&,int);
    void gcsl_neumann_hy(slice&,int);
    void gcsl_neumann_x(slice&,int);
    void gcsl_neumann_int(sliceint&,int);
    void gcsl_neumann_V_int(int*,int);
	void gcsl_noslip(slice&,int);
    void gcsl_outflow(slice&,int);
    void gcsl_potentialbc(slice&,int);
    void gcsl_sommerfeld(slice&,int);

    enum bc_labels { NONE=0, DIRICHLET_ORTH=1, NEUMANN=4, NOSLIP=5, OUTFLOWBC=6, SOMMERFELD=7, POTENTIAL=8,
         DIRICHLET_ORTH_REFLECT=11, DIRICHLET_PARA_REFLECT=12, NEUMANN_X=14, NEUMANN_HX=41, NEUMANN_HY=42,
         HEATBC=61};

    enum dir_labels { X_NEG=1, X_POS=4, Y_NEG=3, Y_POS=2, Z_NEG=5, Z_POS=6 };

    enum gbc_labels { INFLOW=1, OUTFLOW, SYMMETRY, WAVEGEN=6, NUMBEACH, WALL=21};

    // 3D BC decision
    void gcdistro1(field&, int, int, int, int, int, int);
    void gcdistro2(field&, int, int, int, int, int, int);
    void gcdistro3(field&, int, int, int, int, int, int);
    void gcdistro4(field&, int, int, int, int, int, int);

    void gcdistro(field&, int, int, int, bc_labels, int);

    bc_labels gceval1(lexer*,int,int,int);
    bc_labels gceval2(lexer*,int,int,int);
    bc_labels gceval3(lexer*,int,int,int);
    bc_labels gceval4(lexer*,int,int,int);

    // 2D BC decision
    void gcsldistro1(lexer*, slice&, int, int, int, int, int);
    void gcsldistro2(lexer*, slice&, int, int, int, int, int);
    void gcsldistro4(lexer*, slice&, int, int, int, int, int);
    void gcsldistro4a(lexer*, slice&, int, int, int, int, int);

    void gcsldistro1int(sliceint&, int, int, int);
    void gcsldistro2int(sliceint&, int, int, int);
    void gcsldistro4int(sliceint&, int, int, int);
    void gcsldistro4Vint(lexer*, int*, int, int, int);

    void gcsldistro(lexer*, slice&, int, int, bc_labels, int);

    bc_labels gcsleval1(int,int,int);
    bc_labels gcsleval2(int,int,int);
    bc_labels gcsleval4(int,int,int);

    template<typename FlagT, typename GcdfT>
    void gcdf_update_impl(lexer *p, FlagT &flagsf, GcdfT &gcdf, int &gcdf_count);

    MPI_Comm cart_comm = MPI_COMM_NULL;
    int neighbors[6] = {MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL,
                        MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL};
    bool do_comms = true;
    int ndims;

    const int margin;
    int paramargin;
    int m,q,qq;
    double gamma;
    int orderdir;
    int count;
    double starttime,endtime;

    bc_labels bc_label; 
    bc_labels gcval_topodist;
    bool gclabel_outflow;
    bc_labels gclabel_u, gclabel_v, gclabel_w;
    bc_labels gclabel_u_orth,gclabel_v_orth,gclabel_w_orth;
    bc_labels gclabel_u_in,gclabel_v_in,gclabel_w_in,gclabel_press_in,gclabel_lsm_in;
    bc_labels gclabel_u_out, gclabel_v_out, gclabel_w_out;
    bool awa_label,pressout_label,pressin_label;

    // PARALLEL
    double *send1,*send2,*send3,*send4,*send5,*send6;
    double *recv1,*recv2,*recv3,*recv4,*recv5,*recv6;
    int *isend1,*isend2,*isend3,*isend4,*isend5,*isend6;
    int *irecv1,*irecv2,*irecv3,*irecv4,*irecv5,*irecv6;
    double recvsum,recvmin,recvmax;
    int recvisum,recvimin,recvimax;

    int nb0[6],stag[6],rtag[6];

    MPI_Request sreq[6],rreq[6];
    MPI_Status status;

    MPI_Request sreq1,sreq2,sreq3,sreq4,sreq5,sreq6;
    MPI_Request rreq1,rreq2,rreq3,rreq4,rreq5,rreq6;

    double v1,v2,v3,v4;
    double wa,wb;
    double x1,x2;
    double value;

    lexer *p;
    fdm *a;
    fdm2D *b;
    fdm_fnpf *c;
    fdm_nhf *d;
};
#endif
