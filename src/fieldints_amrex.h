/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef FIELDINTS_AMReX_H_
#define FIELDINTS_AMReX_H_

#include"fieldint_amrex.h"

class fieldint1 : public fieldint_amrex
{
public:
    fieldint1(lexer *p) : fieldint_amrex(p) {};
    virtual ~fieldint1() = default;
};

class fieldint2 : public fieldint_amrex
{
public:
    fieldint2(lexer *p) : fieldint_amrex(p) {};
    virtual ~fieldint2() = default;
};

class fieldint3 : public fieldint_amrex
{
public:
    fieldint3(lexer *p) : fieldint_amrex(p) {};
    virtual ~fieldint3() = default;
};

class fieldint4 : public fieldint_amrex
{
public:
    fieldint4(lexer *p) : fieldint_amrex(p) {};
    virtual ~fieldint4() = default;
};

/*!
 * @brief Integer counterpart of field7: sigma-grid vertical-node layout.
 *
 * NODE_Z storage -- z-nodal BoxArray, one plane more than there are cells --
 * so this addresses the same index space as fields_amrex.h's field7 and as the
 * FIJK family does on the flat side.
 *
 * Unlike field7 this carries NO boundary-condition machinery: no
 * FillDomainBoundary, no UpdateBCRecs, no BcDecision. The internal exchange is
 * all these ints need, and fieldint_amrex::FillBoundary OverrideSyncs the
 * shared z-seam plane before filling ghosts.
 */
class fieldint7 : public fieldint_amrex
{
public:
    fieldint7(lexer *p) : fieldint_amrex(p, DataLocation::NODE_Z) {};
    virtual ~fieldint7() = default;
};

#endif
#endif
