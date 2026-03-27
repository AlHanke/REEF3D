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
#ifndef FIELDINT_AMREX_H_
#define FIELDINT_AMREX_H_

#include "fieldint.h"
#include "lexer.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_Vector.H>

class fieldint_amrex : public fieldint
{
public:
    virtual ~fieldint_amrex() = default;

    inline int& operator()(int ii, int jj, int kk) noexcept override final
    {
        return (mf[p->level][*(p->amr_cell_mfi)].array()(amrex::IntVect(AMREX_D_DECL(ii, jj, kk)) + amrex::IntVect(amrex::IntVect(p->amr_tile_lo)), 0));
    }

    inline const int& operator()(int ii, int jj, int kk) const noexcept override final
    {
        return (mf[p->level][*(p->amr_cell_mfi)].const_array()(amrex::IntVect(AMREX_D_DECL(ii, jj, kk)) + amrex::IntVect(amrex::IntVect(p->amr_tile_lo)), 0));
    }

    inline int& operator()(const amrex::IntVect& iv, int comp = 0) noexcept override final
    {
        return (mf[p->level][*(p->amr_cell_mfi)].array()(iv, comp));
    }

    inline const int& operator()(const amrex::IntVect& iv, int comp = 0) const noexcept override final
    {
        return (mf[p->level][*(p->amr_cell_mfi)].const_array()(iv, comp));
    }

    void setVal(int val, bool includeGhost = false) override final;

    void FillBoundary() override final;

    inline amrex::iMultiFab& GetMultiFab() {return mf[p->level];};
    inline const amrex::iMultiFab& GetMultiFab() const {return mf[p->level];};
    inline amrex::iMultiFab& GetMultiFab(int level) {return mf[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return mf[p->level];};

protected:
    fieldint_amrex(lexer* p);

    lexer *p = nullptr;
    amrex::Vector<amrex::iMultiFab> mf = {};
};

#endif
#endif
