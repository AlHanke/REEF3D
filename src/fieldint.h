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

#ifndef FIELDINT_H_
#define FIELDINT_H_

#if not USE_AMREX
#include "lexer.h"
#include <vector>
#endif

#include"field_base.h"
class fieldint : public field_base<int>
{
public:
    virtual ~fieldint() = default;
    #if not USE_AMREX
    fieldint(lexer* p) : p(p), data((p->imax-p->imin)*(p->jmax-p->jmin)*(p->kmax-p->kmin), 0.0) {};
    int& operator()(int ii, int jj, int kk) override {return data[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];};
    void setVal(int val, bool includeGhost = false ) override {int i,j,k;if(includeGhost){MALOOP{operator()(i,j,k) = val;}}else{LOOP{operator()(i,j,k) = val;}}};
private:
    lexer* p;
    std::vector<int> data;
    #endif
};

#endif
