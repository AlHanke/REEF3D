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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef SEDPART_H_
#define SEDPART_H_

#include "sediment.h"
#include "particle_func.h"
#include"increment.h"

#include "particles_obj.h"
#include "field4.h"


class lexer;
class fdm;
class ghostcell;
class ioflow;
class solver;
class reinitopo;
class fdm2D;
class slice;
class ofstrem;
class vrans;
class turbulence;
class bedshear;

/// @brief Class handling sediment on particle basis
/// This class used particles on a lagrangien framework and a VRANS sediment domain to simulate the influence of flow on the sediment
class sedpart : public sediment, private particle_func, private increment
{
public:
    sedpart(lexer*,ghostcell*,turbulence*);
    virtual ~sedpart();

    void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*) override;
    void ini_cfd(lexer*,fdm*,ghostcell*) override;
    void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*) override;
    void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*) override;
    
    void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&) override;
    void ini_sflow(lexer*, fdm2D*, ghostcell*) override;
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*) override;
    
    // ---
    void erode(lexer*,fdm*,ghostcell*);
    // ---
	
    void relax(lexer*,ghostcell*) override;
	double bedshear_point(lexer*,fdm*,ghostcell*) override;
    
    double qbeval(int,int) override;
    void qbeget(int,int,double) override;
    
    double bedzhval(int,int) override;
    
    void ctimesave(lexer*, fdm*) override;
    
    void print_2D_bedload(lexer*, ghostcell*,ofstream&) override;
    void print_3D_bedload(lexer*, ghostcell*,ofstream&) override;
	void name_pvtu_bedload(lexer*, ghostcell*,ofstream&) override;
    void name_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtp_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override;
    
	void print_2D_bedshear(lexer*, ghostcell*,ofstream&) override;
    void print_3D_bedshear(lexer*, ghostcell*,ofstream&) override;
	void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&) override;
    void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtp_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override;
    
    void print_2D_parameter1(lexer*, ghostcell*,ofstream&) override;
    void print_3D_parameter1(lexer*, ghostcell*,ofstream&) override;
	void name_pvtu_parameter1(lexer*, ghostcell*,ofstream&) override;
    void name_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtp_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override;
    
    void print_2D_parameter2(lexer*, ghostcell*,ofstream&) override;
    void print_3D_parameter2(lexer*, ghostcell*,ofstream&) override;
	void name_pvtu_parameter2(lexer*, ghostcell*,ofstream&) override;
    void name_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtp_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override;
    void offset_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override;

protected:

private:
    int maxparticle;
    int gparticle_active;
    int gremoved;
    int gxchange;

    particles_obj PP;
    
    #define PARTLOOP for(int n=0;n<PP.loopindex;++n)
    void seed_ini(lexer*,fdm*,ghostcell*);
    void seed(lexer*,fdm*,ghostcell*);
    void posseed_box(lexer*,fdm*,ghostcell*);
    void posseed_topo(lexer*,fdm*,ghostcell*);
    void posseed_suspended(lexer*,fdm*,ghostcell*);
    void allocate(lexer*,fdm*,ghostcell*);

    field4 active_box;
	field4 active_topo;
    int ppcell;
    int partnum;
    int gpartnum;
    const int irand;
	const double drand;

    // PRINT
	void print_particles(lexer*,fdm*,ghostcell*);
	void print_vtu(lexer*,fdm*,ghostcell*);
	
	void pvtu_pos(lexer*,fdm*,ghostcell*);
    void header_pos(lexer*,fdm*,ghostcell*);
    void piecename_pos(lexer*,fdm*,ghostcell*,int);
	
	char name[100];
    char pname[100];
    double printtime;
    int printcount;
    int num;
    int* changed;
    double* change;

    vrans* pvrans;
    bedshear* pbedshear;
};

#endif