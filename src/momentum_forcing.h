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

#ifndef MOMENTUM_FORCING_H_
#define MOMENTUM_FORCING_H_

#include"increment.h"
#include<vector>

class lexer;
class fdm;
class ghostcell;
class field;
class sixdof;
class fsi;

using namespace std;

class momentum_forcing : public increment
{
public:
    momentum_forcing() = default;
    virtual ~momentum_forcing() = default;
    void momentum_forcing_start(fdm*,lexer*,ghostcell*, sixdof*, fsi*,
                                field&,field&,field&,field&,field&,field&,int,double,bool);

private:
    static constexpr int gcval_u=10, gcval_v=11, gcval_w=12;
};
#endif
