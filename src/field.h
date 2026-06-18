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

#ifndef FIELD_H_
#define FIELD_H_

#include "field_base.h"

#if USE_AMREX
    namespace amrex
    {
        class MultiFab;
    }
    #include "AMReX_Array4.H"
    struct LocalArr4Const {
        amrex::Array4<amrex::Real const> arr;
        int ox, oy, oz; // = tile_lo.{x,y,z}
        AMREX_FORCE_INLINE LocalArr4Const(amrex::Array4<amrex::Real const> a, int x, int y, int z) noexcept
            : arr(a), ox(x), oy(y), oz(z) {}
        AMREX_FORCE_INLINE double operator()(int i, int j, int k) const noexcept {
            return arr(i + ox, j + oy, k + oz);
        }
    };
    struct LocalArr4 {
        amrex::Array4<amrex::Real> arr;
        int ox, oy, oz; // = tile_lo.{x,y,z}
        AMREX_FORCE_INLINE LocalArr4(amrex::Array4<amrex::Real> a, int x, int y, int z) noexcept
            : arr(a), ox(x), oy(y), oz(z) {}
        AMREX_FORCE_INLINE double& operator()(int i, int j, int k) const noexcept {
            return arr(i + ox, j + oy, k + oz);
        }
    };
    struct LocalArr4IntConst {
        amrex::Array4<int const> arr;
        int ox, oy, oz;
        AMREX_FORCE_INLINE LocalArr4IntConst(amrex::Array4<int const> a, int x, int y, int z) noexcept
            : arr(a), ox(x), oy(y), oz(z) {}
        AMREX_FORCE_INLINE int operator()(int i, int j, int k) const noexcept {
            return arr(i + ox, j + oy, k + oz);
        }
    };
#endif

class field : public field_base<double>
{
public:
#if USE_AMREX
    field() = default;
#else
    field(lexer* p) : field_base<double>(p) {}
#endif
    virtual ~field() = default;

#if USE_AMREX
    virtual void FillDomainBoundary(int gcv) = 0;
    virtual void FillDomainBoundaryValue(double value, int dir, bool high) = 0;
    virtual amrex::MultiFab& GetMultiFab() = 0;
    virtual const amrex::MultiFab& GetMultiFab() const = 0;
    virtual amrex::MultiFab& GetMultiFab(int) = 0;
    virtual const amrex::MultiFab& GetMultiFab(int) const = 0;
#endif

#if not USE_AMREX
    void CopyFrom(const field& src) {V = src.V; cache_addressing();};
#endif
};

#endif
