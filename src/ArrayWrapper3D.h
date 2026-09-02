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

/*!
    * @brief Values for the data_location constructor argument.
    *
    * 1/2/3 mirror DataLocation::FACE_X/Y/Z and select the
    * staggered direction for the boundary shift; 4 is cell-centred, the
    * default. 7 is the sigma-grid vertical-node layout: vertical stride
    * p->kmaxF instead of p->kmax, one plane more than there are cells, which
    * reproduces the FIJK addressing in iterators3D.h. It is the integer
    * counterpart of field7 and the "7" in flag7 / gcx_flagx7 / gcx_parax7co.
    *
    * The extent is a constructor argument rather than a distinct type on
    * purpose: this class is final and its accessors are hot (see
    * cache_addressing on why the strides are long), so a subclass would cost
    * a vtable on every flag test for no expressive gain.
*/

class ArrayWrapper3D final
{
public:
    explicit ArrayWrapper3D(lexer *pp, DataLocation _data_location = DataLocation::CELL_CENTERED);
    #if USE_AMREX
    /// View constructor: non-owning view into component @p comp of @p shared.
    /// The caller must ensure @p shared outlives this object.
    explicit ArrayWrapper3D(lexer *pp, amrex::Vector<amrex::iMultiFab> *shared, int comp);
    #endif

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

    #if not USE_AMREX
    /// Lightweight, capture-by-value view — mirrors field_base::View/ConstView
    /// (see field_base.h for the rationale: base pointer + strides read once
    /// instead of through operator()'s `this->m_base/m_js/m_ks` on every call,
    /// so a hoisted flag check survives opaque calls inside the loop body).
    struct ConstView
    {
        const int* __restrict base;
        long js, ks;
        inline const int& operator()(int i, int j, int k) const noexcept
        { return base[i*js + j*ks + k]; }
    };

    inline ConstView const_view() const noexcept { return ConstView{m_base, m_js, m_ks}; }
    #endif

    #if USE_AMREX
    void fillBoundary();
    /// Batch FillBoundary: fills @p ncomp components starting at @p scomp of
    /// @p shared in a single MPI collective instead of one per component.
    static void FillBoundaryBatch(lexer* p, amrex::Vector<amrex::iMultiFab>& shared,
                                  int scomp, int ncomp);
    void fillHigherLevels();

    inline amrex::iMultiFab& GetMultiFab();
    inline const amrex::iMultiFab& GetMultiFab() const;

    inline amrex::iMultiFab& GetMultiFab(int level) {return m_shared ? (*m_shared)[level] : data[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return m_shared ? (*m_shared)[level] : data[level];};
    #endif

private:
    /// Vertical extent for this wrapper's data location: p->kmaxF for NODE_Z,
    /// p->kmax otherwise. Resolved on demand rather than cached at construction
    /// — grid::assign_margin sets kmax and kmaxF together and both are final by
    /// the time resize() runs, but not necessarily by the time the ctor does
    /// (lexer builds its own wrappers in its initialiser list).
    int kz() const noexcept;

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

    amrex::Vector<amrex::iMultiFab>* m_shared = nullptr; ///< non-owning ptr (view mode only)
    int m_comp = 0;                                       ///< component in m_shared

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
