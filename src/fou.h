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

#ifndef FOU_H_
#define FOU_H_

#include"increment.h"
#include"convection.h"

class flux;

class fou final : public convection, public increment
{
public:
    fou(lexer*);
    virtual ~fou() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    template<typename GenericField>
    inline double aij(lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&, double*, double*, double*);

    flux *pflux;
};

#endif
