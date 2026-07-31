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

#ifndef KEPSILON_BC_H_
#define KEPSILON_BC_H_

#include "increment.h"
#include "roughness.h"

class lexer;
class fdm;
class field;

class kepsilon_bc : public roughness
{
public:
    kepsilon_bc(lexer*);
    virtual ~kepsilon_bc() = default;
    void bckeps_start(fdm*,lexer*,field&,field&,int);
    void bckin_matrix(lexer*,fdm*,field&);
    void bcepsilon_matrix(lexer*,fdm*,field&);

private:
    void wall_law_kin(lexer*,fdm*,field&,field&,int,int,int,int,int,int);
    void wall_law_eps(lexer*,fdm*,field&,field&,int,int,int,int,int,int);
};
#endif
