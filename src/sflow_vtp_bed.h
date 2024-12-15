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

#ifndef SFLOW_VTP_BED_H_
#define SFLOW_VTP_BED_H_

#include"increment.h"
#include"vtp3D.h"
#include<vector>

class lexer;
class fdm2D;
class ghostcell;
class sediment;

class sflow_vtp_bed : virtual public increment, private vtp3D
{
public:
    sflow_vtp_bed(lexer*,fdm2D*);
    virtual ~sflow_vtp_bed();
    
    void start(lexer*,fdm2D*,ghostcell*,sediment*);
    
private:
    void print2D(lexer*,fdm2D*,ghostcell*,sediment*);
    void pvtp(lexer*,fdm2D*,ghostcell*,sediment*);
    void name_iter(lexer*,fdm2D*,ghostcell*);
    void piecename(lexer*,fdm2D*,ghostcell*,int);
    
    std::vector<char> buffer;
    int m;
    
    char name[200],pname[200];
    int n,iin,offset[200];
    float ffn;
    
    int printbedcount;
    double printbedtime;
};

#endif
