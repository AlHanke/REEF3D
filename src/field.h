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

#ifndef FIELD_H_
#define FIELD_H_

#include "field_base.h"

#if not USE_AMREX
#include "lexer.h"
#include <vector>
#include <algorithm>
#endif

class field : public field_base<double>
{
public:
    virtual ~field() = default;
    #if USE_AMREX
    virtual void FillDomainBoundary(int gcv) = 0;
    virtual void FillDomainBoundaryValue(double value, int dir, bool high) = 0;
    #else
    field(lexer* p) : p(p), data((p->imax-p->imin)*(p->jmax-p->jmin)*(p->kmax-p->kmin), 0.0) {};
    double& operator()(int ii, int jj, int kk) override final {return data[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin];};
    void setVal(double val, bool includeGhost = false ) override final {int i,j,k;if(includeGhost){std::fill(data.begin(),data.end(),val);}else{LOOP{operator()(i,j,k) = val;}}};
private:
    lexer* p;
    std::vector<double> data;
    #endif
};

#endif
