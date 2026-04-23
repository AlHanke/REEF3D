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

#ifndef VANALBADA_H_
#define VANALBADA_H_

#include"fluxlim.h"
#include"increment.h"

class vanalbada final : public fluxlim, public increment
{
public:
    vanalbada(lexer*) {};
    virtual ~vanalbada() = default;

protected:
    inline double phi_impl(double vn1, double vn2, double vq1, double vq2, double) override final
    {
        double denom = vq1 - vq2;
        double r = (vn1 - vn2) / (fabs(denom)>1.0e-10 ? denom : 1.0e20);
        return (r*r + r) / (r*r + 1.0);
    };
};

#endif

