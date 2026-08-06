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

#ifndef ARRAYWRAPPER3D_H_
#define ARRAYWRAPPER3D_H_

#include "definitions.h"

#if USE_AMREX
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_Vector.H>
#else
#include <vector>
#endif

class lexer;

class ArrayWrapper3D final
{
public:
    explicit ArrayWrapper3D(lexer *pp, DataLocation _data_location = DataLocation::CELL_CENTERED);

    ArrayWrapper3D(const ArrayWrapper3D&)            = delete;
    ArrayWrapper3D &operator=(const ArrayWrapper3D&) = delete;
    ArrayWrapper3D(ArrayWrapper3D&&)                 = delete;
    ArrayWrapper3D &operator=(ArrayWrapper3D&&)      = delete;

    ~ArrayWrapper3D();

    void resize(int default_value = 0);

    void setVal(int val, bool includeGhost = false);

    inline int &operator()(int i, int j, int k) noexcept;
    inline const int &operator()(int i, int j, int k) const noexcept;

    inline int &operator[](int index) noexcept;
    inline const int &operator[](int index) const noexcept;

    operator int *();

    #if USE_AMREX
    void fillBoundary();
    void fillHigherLevels();
    amrex::iMultiFab& GetMultiFab();
    const amrex::iMultiFab& GetMultiFab() const;
    inline amrex::iMultiFab& GetMultiFab(int level) {return data[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return data[level];};
    #endif

private:
    #if not USE_AMREX
    /// Recomputes the cached addressing below from @p data and the lexer grid
    /// metrics. Must be called after anything that resizes @p data.
    void cache_addressing() noexcept;
    #endif

    lexer *p = nullptr;

    DataLocation data_location = DataLocation::CELL_CENTERED;

    #if USE_AMREX
    /// Memoizes the Array4 and tile origin for the fab the lexer currently points
    /// at. Const so both the const and non-const accessors can share one cache;
    /// the state below is mutable because filling it is logically const.
    AMREX_FORCE_INLINE void refresh_cache_if_needed() const noexcept;

    amrex::Vector<amrex::iMultiFab> data;

    mutable amrex::Array4<int> m_cached_arr4 = {};
    mutable int m_cached_ox      = 0;
    mutable int m_cached_oy      = 0;
    mutable int m_cached_oz      = 0;
    mutable int m_cached_mfi_idx = -1;
    mutable int m_cached_level   = -1;
    mutable int m_cached_til_idx = -1;
    #else
    std::vector<int> data;

    int* m_base = nullptr; ///< origin-folded base: m_base[i*m_js + j*m_ks + k]
    long m_js   = 0;       ///< i-stride (jmax*kmax); long on purpose, see cache_addressing
    long m_ks   = 0;       ///< j-stride (kmax);      long on purpose, see cache_addressing
    #endif
};

#endif
