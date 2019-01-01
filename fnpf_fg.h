/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

class lexer;
class fdm;
class ghostcell;
class solver;
class convection;
class ioflow;
class reini;
class onephase;

using namespace std;

#ifndef FNPF_FG_H_
#define FNPF_FG_H_

class fnpf_fg
{
public:
	virtual void start(lexer*, fdm*, ghostcell*, solver*, convection*, ioflow*, reini*,onephase*)=0;
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*, reini*, onephase*)=0;
    virtual void inidisc(lexer*, fdm*, ghostcell*)=0;

};

#endif
