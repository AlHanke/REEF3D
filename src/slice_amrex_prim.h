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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef SLICE_AMREX_PRIM_H_
#define SLICE_AMREX_PRIM_H_

class slice_amrex_prim
{
public:
    slice_amrex_prim() = default;
    virtual ~slice_amrex_prim() = default;

    virtual void regrid() = 0;
};

#endif
#endif
