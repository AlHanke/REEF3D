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

#ifndef NHFLOW_VTP_FSF_H_
#define NHFLOW_VTP_FSF_H_

#include"increment.h"
#include"vtp3D.h"
#include<vector>

class lexer;
class fdm_nhf;
class ghostcell;
class ioflow;
class sediment;

using namespace std;

class nhflow_vtp_fsf : virtual public increment, private vtp3D
{
public:
	nhflow_vtp_fsf(lexer*,fdm_nhf*,ghostcell*);
	virtual ~nhflow_vtp_fsf() = default;
	
    void start(lexer*,fdm_nhf*,ghostcell*,sediment*);
    void print2D(lexer*,fdm_nhf*,ghostcell*,sediment*);
    void preproc(lexer*,fdm_nhf*,ghostcell*);
	
private:
	void pvtu(lexer*,fdm_nhf*,ghostcell*,sediment*);
	void name_iter(lexer*,fdm_nhf*,ghostcell*);
    void piecename(lexer*,fdm_nhf*,ghostcell*,int);
	
    std::vector<char> buffer;
    int m = 0;
	
	char name[200],pname[200];
    int n,iin,offset[200];
    float ffn;
    
	double xs_local,ys_local,zs_local,xe_local,ye_local,ze_local;
	double xs_global,ys_global,zs_global,xe_global,ye_global,ze_global;
    
    int gcval_eta, gcval_fifsf;
    int printcount;
    int jj;
    
    int *wetmax;
    
    ofstream pvdout;
	

};

#endif
