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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#ifndef FIELD_BASE_H_
#define FIELD_BASE_H_

#if not USE_AMREX
    #include "lexer.h"
    #include <vector>
#endif

template<typename T>
class field_base
{
public:
#if USE_AMREX
    field_base() = default;
#else
    field_base(lexer* p) :
        imin(p->imin), imax(p->imax),
        jmin(p->jmin),
        kmin(p->kmin), kmax(p->kmax),
        jkmax(p->jmax*p->kmax),
        p(p),
        V(imax*jkmax, T{})
    {};
#endif
    field_base(const field_base&) = delete;
    field_base& operator=(const field_base&) = delete;
    field_base(field_base&&) = delete;
    field_base& operator=(field_base&&) = delete;
    virtual ~field_base() = default;

#if USE_AMREX
    virtual T& operator()(int, int, int) noexcept = 0;
#else
    inline T& operator()(int ii, int jj, int kk) noexcept {return V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin];};
#endif

#if USE_AMREX
    virtual void FillBoundary() = 0;
#endif

#if not USE_AMREX
private:
    const int imin,imax,jmin,kmin,kmax,jkmax;
    lexer *p;

protected:
    std::vector<T> V;
#endif
};

#endif
