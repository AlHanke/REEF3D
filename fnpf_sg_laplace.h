/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

class lexer;
class fdm_fnpf;
class field;
class ghostcell;
class solver;
class vec;
class fnpf_sg_fsfbc;

#ifndef FNPF_SG_LAPLACE_H_
#define FNPF_SG_LAPLACE_H_

using namespace std;

class fnpf_sg_laplace
{
public:

    virtual void start(lexer *,fdm_fnpf*,ghostcell*,solver*,fnpf_sg_fsfbc*,double*)=0;

};

#endif
