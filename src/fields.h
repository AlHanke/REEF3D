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

#if not USE_AMREX
#ifndef FIELDS_H_
#define FIELDS_H_

#include "field.h"

class field1 : public field
{
public:
    field1(lexer* p) : field(p) {};
    virtual ~field1() = default;
};

class field2 : public field
{
public:
    field2(lexer* p) : field(p) {};
    virtual ~field2() = default;
};

class field3 : public field
{
public:
    field3(lexer* p) : field(p) {};
    virtual ~field3() = default;
};

class field4 : public field
{
public:
    field4(lexer* p) : field(p) {};
    virtual ~field4() = default;
};

/*!
 * @brief Sigma-grid vertical-node field: n*m*(o+1) planes, vertical stride p->kmaxF.
 *
 * Same (i,j,k) accessor as field1-4 — only the vertical extent differs. The
 * vertical index runs 0..KMAX_LOOP+1 (FKLOOP) rather than 0..KMAX_LOOP (KLOOP),
 * i.e. one plane per cell plus a closing top plane, which is what the sigma-grid
 * solvers (FNPF, NHFLOW) need. The flat-array equivalent is the FIJK family in
 * iterators3D.h; field7 exists so those call sites can drop the flat index and
 * address by (i,j,k) like every other field.
 */
class field7 : public field
{
public:
    /*!
     * The slack plane mirrors the imax*jmax*(kmax+2) allocation in
     * driver_makegrid_sigma.cpp against a kmaxF = kmax+1 stride. FIJKp3 at the
     * top node of the last column (fnpf_fsf_update.cpp, inside FFILOOP4) lands
     * on the very last in-stride slot, so the spare plane is the only headroom
     * the forward-stencil macros have. Drop it only after auditing FIJKp3/p4.
     */
    field7(lexer* p) : field(p, p->kmaxF,
                             static_cast<std::size_t>(p->imax)*p->jmax) {};
    virtual ~field7() = default;

    /*!
     * FBASELOOP rather than LOOP: LOOP stops at KMAX_LOOP and would leave the
     * top node plane untouched, and its PCHECK reads the IJK-strided flag4.
     * Shadows the non-virtual field_base::setVal — field7 is meant to be held
     * by concrete type, never through a field& (CopyFrom would size-mismatch).
     */
    void setVal(double val, bool includeGhost = false)
    {
        if(includeGhost)
        {
            std::fill(V.begin(), V.end(), val);
            return;
        }

        int i,j,k;
        FBASELOOP
        operator()(i,j,k) = val;
    };
};

#endif
#endif
