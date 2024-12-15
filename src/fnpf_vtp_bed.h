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

#ifndef FNPF_VTP_BED_H_
#define FNPF_VTP_BED_H_

#include"increment.h"
#include"vtp3D.h"
#include<vector>

class lexer;
class fdm_fnpf;
class ghostcell;
class ioflow;

class fnpf_vtp_bed : public increment, private vtp3D
{
public:
	fnpf_vtp_bed(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_vtp_bed();
	
    void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*);
	
private:
	void print2D(lexer*,fdm_fnpf*,ghostcell*);
	void pvtp(lexer*,fdm_fnpf*,ghostcell*);
	void name_iter(lexer*,fdm_fnpf*,ghostcell*);
    void piecename(lexer*,fdm_fnpf*,ghostcell*,int);

    std::vector<char> buffer;
    int m = 0;
	
	char name[200],pname[200];
    int n,iin,offset[200];
    float ffn;

    int printcount;

};

#endif
