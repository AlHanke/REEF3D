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

#ifndef ARRAYWRAPPER2D_H_
#define ARRAYWRAPPER2D_H_

// =====================================================================
// ArrayWrapper2D — 2D (column) integer array for flagslice1/2/4.
//
// Relationship to ArrayWrapper3D mirrors slice_amrex : field_amrex.
//   - data lives on z-slabs (makeSlab(box3d, 2, 0)) of the 3D BoxArray
//   - ghost cells in x/y only: IntVect(margin, margin, 0)
//   - view/unique/owner split, because the slab layout is OVERLAPPING
//     (every 3D box in a column flattens onto the same (i,j)) and an
//     overlapping MultiFab can never be a FillBoundary/ParallelCopy source
//
// This header is included into lexer.h BEFORE the lexer class is formed, so
// it must not reference any lexer member. It therefore only DECLARES the
// methods that touch `p->...`:
//   - inline hot-path methods (operator(), operator[], refresh caches,
//     GetMultiFab()) are defined in ArrayWrapper2D_imp.h, included
//     at the bottom of lexer.h once lexer is complete (same trick as
//     ArrayWrapper3D_imp.h).
//   - everything else is defined in ArrayWrapper2D.cpp.
//
// Two things differ from slice_amrex, both because the payload is int:
//
//   1. No FillPatch / no BCRec / no MyExtBCFillSlice. amrex::Interpolater and
//      FillPatchTwoLevels are amrex::Real-only, so coarse->fine goes through
//      the same manual injection ArrayWrapper3D::fillHigherLevels() uses,
//      with the z ratio pinned to 1. Domain-boundary flags stay owned by the
//      gcsl_* / mgcslice* code, exactly as flag1-4's ghosts are owned by
//      ghostcell::start1-4 rather than by ArrayWrapper3D itself.
//
//   2. The owner mask is not ksurf-based. A flag is a column property with no
//      surface to straddle, so ownership goes to the box touching the domain
//      z-lo face. Uniqueness is free: two boxes at one level sharing z-lo and
//      overlapping in (i,j) would overlap in 3D, which a BoxArray forbids. On
//      fine levels that do not reach the bottom a column has NO owner — that is
//      a hole, filled by injection from lev-1, same as slice_amrex::fillHoles().
// =====================================================================

#if USE_AMREX
#include "slice_amrex_prim.h"
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#else
#include <vector>
#endif

class lexer;

class ArrayWrapper2D final
#if USE_AMREX
    : public slice_amrex_prim
#endif
{
public:
    explicit ArrayWrapper2D(lexer *pp);

    ArrayWrapper2D(const ArrayWrapper2D&)            = delete;
    ArrayWrapper2D &operator=(const ArrayWrapper2D&) = delete;
    ArrayWrapper2D(ArrayWrapper2D&&)                 = delete;
    ArrayWrapper2D &operator=(ArrayWrapper2D&&)      = delete;

    ~ArrayWrapper2D();

    // Deferred define, like ArrayWrapper3D::resize(). NOT register_imf():
    // that registry re-defines its entries with the 3D amrex_box_array on
    // regrid, which would silently un-flatten the slabs. Slabs go through
    // slice_registry and rebuild themselves in regrid() below.
    // -10 is a sentinel value for "uninitialized" (the default)
    void resize(int default_value = -10);

    void setVal(int val, bool includeGhost = false);

    inline int &operator()(int ii, int jj) noexcept;
    inline const int &operator()(int ii, int jj) const noexcept;

    // IJ-style flat access, level 0 only — the same escape hatch (and the same
    // restriction) as ArrayWrapper3D::operator[]. IJ is
    // (i-imin)*jmax + (j-jmin), so decode against jmax.
    inline int &operator[](int index) noexcept;
    inline const int &operator[](int index) const noexcept;

    operator int *();

    #if USE_AMREX
    void fillBoundary();
    void fillHigherLevels();
    void regrid() override final;

    inline amrex::iMultiFab& GetMultiFab();
    inline const amrex::iMultiFab& GetMultiFab() const;
    inline amrex::iMultiFab& GetMultiFab(int level) noexcept { return m_view[level]; };
    inline const amrex::iMultiFab& GetMultiFab(int level) const noexcept { return m_view[level]; };
    #endif

private:
    #if USE_AMREX
    void define_level(int lev, bool seed_from_coarse);
    void buildOwnerMask(int lev);
    void makeUnique(int lev);
    void fillHoles(int lev);

    // Piecewise-constant injection of coarse (fine_lev-1) unique data onto
    // fine_mf's valid + x/y ghost cells that lie inside the fine domain. Shared
    // by fillHigherLevels() (fine_mf = m_view) and define_level()'s coarse seed
    // (fine_mf = the fresh unique MultiFab). z ratio is pinned to 1.
    void inject_from_coarse(amrex::iMultiFab& fine_mf,
                            const amrex::iMultiFab& crse_unique,
                            int fine_lev);

    AMREX_FORCE_INLINE void refresh_cache_if_needed() noexcept;
    AMREX_FORCE_INLINE void refresh_const_cache_if_needed() const noexcept;

    #else
    /// Recomputes the cached addressing below. Must be called after anything
    /// that resizes @p data. See ArrayWrapper3D::cache_addressing for why the
    /// origin is folded into the pointer and why the stride is a long.
    void cache_addressing() noexcept;
    #endif

    lexer *p = nullptr;

    #if USE_AMREX
    amrex::Vector<amrex::iMultiFab> m_view;     ///< overlapping: what operator() reads
    amrex::Vector<amrex::iMultiFab> m_unique;   ///< removeOverlap partition: comms source
    amrex::Vector<amrex::iMultiFab> m_owner;    ///< 1 on the one view cell per column that counts

    amrex::IntVect m_ghost = amrex::IntVect(AMREX_D_DECL(0,0,0));
    int nlevs     = 0;
    int m_default = -10;   ///< read_grid.cpp's flagslice init value

    amrex::Array4<int> m_cached_arr4 = {};
    int m_cached_ox      = 0;
    int m_cached_oy      = 0;
    int m_cached_mfi_idx = -1;
    int m_cached_level   = -1;
    int m_cached_til_idx = -1;

    mutable amrex::Array4<const int> m_cached_const_arr4 = {};
    mutable int m_cached_const_ox      = 0;
    mutable int m_cached_const_oy      = 0;
    mutable int m_cached_const_mfi_idx = -1;
    mutable int m_cached_const_level   = -1;
    mutable int m_cached_const_til_idx = -1;
    #else
    std::vector<int> data;

    int *m_base = nullptr; ///< origin-folded base: m_base[ii*m_js + jj]
    long m_js   = 0;       ///< i-stride (jmax); long on purpose
    #endif
};

#endif
