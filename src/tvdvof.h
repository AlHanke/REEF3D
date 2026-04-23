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

#ifndef TVDVOF_H_
#define TVDVOF_H_

#include"fluxlim.h"
#include"increment.h"

class tvdvof : public fluxlim, public increment
{

public:

    tvdvof(lexer*) {};
    virtual ~tvdvof() = default;

protected:
    inline double phi_impl(double vn1, double, double, double, double vcell) override final
    {
        const double rp = vcell;
        const double rn = vn1;
        return std::min(std::max(1.0 - pow(std::max(pow(1.0-4.0*rp*(1.0-rp),2.0), 1.0-(1.0-4.0*rn*(1.0-rn))),2.0), 0.0), 1.0);
    };
};

#endif


