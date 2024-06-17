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

#ifndef FNPF_VTS3D_H_
#define FNPF_VTS3D_H_

#include"fnpf_printer.h"
#include"increment.h"

class fdm_fnpf;
class force_ale;
class fnpf_print_wsf;
class fnpf_print_wsf_theory;
class fnpf_print_wsfline;
class fnpf_print_wsfline_y;
class fnpf_vtp_fsf;
class fnpf_vtp_bed;
class fnpf_state;
class fnpf_breaking_log;
class fnpf_print_Hs;
class fnpf_vel_probe;
class fnpf_vel_probe_theory;
class fnpf_runup;
class potentialfile_out;
class ioflow;

using namespace std;

class fnpf_vts3D : public fnpf_printer, public increment
{

public:
	fnpf_vts3D(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_vts3D();
	virtual void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*);
    virtual void print_stop(lexer*,fdm_fnpf*,ghostcell*);
    virtual void print_vtu(lexer*,fdm_fnpf*,ghostcell*);
    
private:
    void pvts(lexer*,ghostcell*);
    void name_iter(lexer*,ghostcell*);
    void name_time(lexer*,ghostcell*);
    void piecename(lexer*,ghostcell*, int);
    void extent(lexer*,int);
    void fextent(lexer*);

    char name[200],pname[200],epsvar[200],pextent[20];
    int iextent[6];
    int *piextent;
    int n,iin,offset[200];
    float ffn;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double *printfsftime_wT;
    int *printfsfiter_wI;
    double phase;
    double zcoor;
    
    int printcount;
    
    fnpf_print_wsf *pwsf;
    fnpf_print_wsf_theory *pwsf_theory;
    fnpf_print_wsfline *pwsfline;
    fnpf_print_wsfline_y *pwsfline_y;
    potentialfile_out *ppotentialfile;
    fnpf_vtp_fsf *pfsf;
    fnpf_vtp_bed *pbed;
    fnpf_state *pstate;
    fnpf_breaking_log *pbreaklog;
	force_ale **pforce_ale;
    fnpf_runup **prunup;
    fnpf_print_Hs *phs;
    fnpf_vel_probe *pvel;
    fnpf_vel_probe_theory *pveltheo;
};

#endif

