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

#ifndef CICSAM_H_
#define CICSAM_H_

#include"convection.h"
#include"increment.h"

class flux;

class cicsam final : public convection, public increment
{
public:
    cicsam(lexer*);
    virtual ~cicsam() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    template<typename GenericField>
    inline double aij(lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&);

    template<typename GenericField>
    double cface(lexer*,fdm*,const GenericField&,int,int,double);

    flux *pflux;
};

#endif
