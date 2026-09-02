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

#ifndef FLUX_H_
#define FLUX_H_

class fdm;
class field;
struct LocalArr4Const;

#if not USE_AMREX
// Needed for field::ConstView below -- a nested type name requires a complete
// type, so the forward declaration above isn't enough on this path.
#include "field.h"
#endif

class flux
{
public:
    virtual ~flux() = default;

    virtual void u_flux(fdm*,int,const field&,double&,double&) = 0;
    virtual void v_flux(fdm*,int,const field&,double&,double&) = 0;
    virtual void w_flux(fdm*,int,const field&,double&,double&) = 0;

    #if USE_AMREX
    virtual void u_flux(fdm*,int,const LocalArr4Const&,double&,double&) = 0;
    virtual void v_flux(fdm*,int,const LocalArr4Const&,double&,double&) = 0;
    virtual void w_flux(fdm*,int,const LocalArr4Const&,double&,double&) = 0;
    #else
    virtual void u_flux(fdm*,int,const field::ConstView&,double&,double&) = 0;
    virtual void v_flux(fdm*,int,const field::ConstView&,double&,double&) = 0;
    virtual void w_flux(fdm*,int,const field::ConstView&,double&,double&) = 0;
    #endif
};

#endif
