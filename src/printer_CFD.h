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

#ifndef printer_CFD_H_
#define printer_CFD_H_

#include"printer.h"
#include"nodefill.h"
#include"field5.h"

#include "vtk3D.h"

#include<sstream>

#include <mpi.h>

class turbulence;
class heat;
class suspended;
class bedload;
class topo;
class print_wsf;
class print_wsf_theory;
class print_wsfline_x;
class print_wsfline_y;
class force;
class force;
class vorticity;
class solver;
class probe_point;
class probe_pressure;
class probe_line;
class bedprobe_point;
class bedprobe_max;
class gage_discharge_x;
class gage_discharge_window_x;
class fsf_vtp;
class topo_vtp;
class cfd_state;
class bedshear_probe;
class bedshear_max;
class sloshing_force;
class print_porous;
class bedprobe_line_x;
class bedprobe_line_y;
class probe_vel;
class probe_vel_theory;
class exportfile;
class flowfile_out;
class print_averaging;

using namespace std;

class printer_CFD : public printer, public nodefill 
{

public:
	printer_CFD(lexer*,fdm*,ghostcell*);
	virtual ~printer_CFD();
	virtual void start(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);
    virtual void print_vtk(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);
    virtual void print_stop(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);

private:
    void print3D(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*);
    void parallelData(fdm*,lexer*,ghostcell*,turbulence*,heat*,data*,concentration*,multiphase*,sediment*);

    /// Compact
    void setupCompactPrint(lexer*,fdm*,ghostcell*);
    void print3Dcompact(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*);

    /// MPI
    void prepMPIprint(lexer*,fdm*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*, sediment*);
    void offsetsComp(lexer*,fdm*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*, sediment*);
    void print3DMPI(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*);
    void writeHeader(lexer*,fdm*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*, sediment*, MPI_File &file, int *piextent);
    void pointData(lexer*,fdm*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*, std::ofstream &result, int *offset);
    void writeData(lexer*,fdm*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*, sediment*, MPI_File &);
    bool MPIinitialPrint=true;
    int *offsets=nullptr;
    MPI_Offset *offsetMPI=nullptr;
    // bool initialMPIprint=true;
    // int offsetMPI[300];

    /// MPIC
    void setupCompactMPIPrint(lexer*,fdm*,ghostcell*);
    void print3DcompactMPI(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*);
    void offsetCMPIPoints(lexer*,int*,int*,int);
    size_t headerSize=0;
    int endIndex=0;
    std::vector<MPI_Offset> offsetCMPI;
    std::vector<int> offsetCMPIitr;
    std::stringstream header;
    std::stringstream footer;
    int kbeginPoint,kendPoint;
    int jbeginPoint,jendPoint;
    int ibeginPoint,iendPoint;
    int *gbeginEndPoint=nullptr;
    int compactMPIPOffset[300];

    double **press;
    double **uvel;
    double **vvel;
    double **wvel;
    double **topo;
    double **phi;
    double **eddyv;
    int **flag;
    int **flag5;
    double *XN;
    double *YN;
    double *ZN;
    int *recvcounts;
    int *displs;
    int *gneibours;
    int *gextent;
    int *globalSendCounts;
    double* pressGlobal;
    double* uvelGlobal;
    double* vvelGlobal;
    double* wvelGlobal;
    double* topoGlobal;
    double* phiGlobal;
    double* eddyvGlobal;
    int* flagGlobal;
    int* flag5Global;
    int localSendCount;
    int cellNum;
    int pointNum;
    int m,compactOffset[300];

    char name[200];
    int n,iin,offset[300];
    float ffn;
    // int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double *printfsftime_wT;
    // int *printfsfiter_wI;
	
	field5 eta;

    print_wsf *pwsf;
	print_wsf_theory *pwsf_theory;
    print_wsfline_x *pwsfline_x;
	print_wsfline_y *pwsfline_y;
    force **pforce;
    int P81;
    vorticity *pvort;
	probe_point *pprobe;
    probe_pressure *ppressprobe;
	probe_line *pline;
	bedprobe_point *pbedpt;
	bedprobe_line_x *pbedlinex;
	bedprobe_line_y *pbedliney;
	bedprobe_max *pbedmax;
	bedshear_probe *pbedshear;
	bedshear_max *pbedshearmax;
	gage_discharge_x *pq;
    gage_discharge_window_x *pqw;
	fsf_vtp *pfsf;
    topo_vtp *ptopo;
	cfd_state *pstate;
    sloshing_force *pslosh;
	print_porous *ppor;
    exportfile *pexport;
    flowfile_out *pflowfile;
    print_averaging *pmean;
    probe_vel *pvel;
    probe_vel_theory *pveltheo;

    vtk3D *outputFormat;
};

#endif

