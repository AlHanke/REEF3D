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
Author: Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef FIELDS_AMReX_H_
#define FIELDS_AMReX_H_

#include "field_amrex.h"

class field1 : public field_amrex
{
public:
    field1 (lexer*);
    virtual ~field1() = default;
    void FillDomainBoundary(int gcv) override;
private:
    amrex_bc_func::Field1BcDecision::Field1Params params;
};

class field2 : public field_amrex
{
public:
    field2 (lexer*);
    virtual ~field2() = default;
    void FillDomainBoundary(int gcv) override;
private:
    amrex_bc_func::Field2BcDecision::Field2Params params;
};

class field3 : public field_amrex
{
public:
    field3 (lexer*);
    virtual ~field3() = default;
    void FillDomainBoundary(int gcv) override;
private:
    amrex_bc_func::Field3BcDecision::Field3Params params;
};

class field4 : public field_amrex
{
public:
    field4 (lexer*);
    virtual ~field4() = default;
    void FillDomainBoundary(int gcv) override;
private:
    amrex_bc_func::Field4BcDecision::Field4Params params;
};

#endif
#endif
