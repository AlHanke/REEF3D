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

#ifndef ARRAYWRAPPER_INT_H_
#define ARRAYWRAPPER_INT_H_

#if USE_AMREX
#include <AMReX_iMultiFab.H>
#include <AMReX_MFIter.H>
#else
#include <vector>
#endif

class lexer;

class ArrayWrapper_int
{
public:
    ArrayWrapper_int(lexer* p);
    virtual ~ArrayWrapper_int() = default;

    void resize(int default_value = 0);

    inline int& operator[] (int index) noexcept;

    #if USE_AMREX
    /// View constructor: non-owning view into component @p comp of @p shared.
    /// The caller must ensure @p shared outlives this object.
    ArrayWrapper_int(lexer* p, amrex::Vector<amrex::iMultiFab>* shared, int comp);

    inline int& operator() (int i, int j, int k) noexcept
    {
        refresh_cache_if_needed();
        return m_cached_arr4(m_cached_ox + i, m_cached_oy + j, m_cached_oz + k, m_comp);
    };

    /// Batch FillBoundary: fills @p ncomp components starting at @p scomp of
    /// @p shared in a single MPI collective instead of one per component.
    static void FillBoundaryBatch(lexer* p, amrex::Vector<amrex::iMultiFab>& shared,
                                  int scomp, int ncomp);
    #endif

    operator int* ();
    void setVal(int val, bool includeGhost = false);

    #if USE_AMREX
    void fillBoundary();
    void fillHigherLevels();
    inline amrex::iMultiFab& GetMultiFab();
    inline const amrex::iMultiFab& GetMultiFab() const;
    inline amrex::iMultiFab& GetMultiFab(int level) {return m_shared ? (*m_shared)[level] : data[level];};
    inline const amrex::iMultiFab& GetMultiFab(int level) const {return m_shared ? (*m_shared)[level] : data[level];};
    #endif

private:
    lexer* p;

    #if USE_AMREX
    void refresh_cache_if_needed() noexcept;

    std::vector<amrex::iMultiFab> data;

    amrex::Vector<amrex::iMultiFab>* m_shared = nullptr; ///< non-owning ptr (view mode only)
    int m_comp = 0;                                       ///< component in m_shared

    amrex::Array4<int> m_cached_arr4 = {};
    int m_cached_ox      = 0;
    int m_cached_oy      = 0;
    int m_cached_oz      = 0;
    int m_cached_mfi_idx = -1;
    int m_cached_level   = -1;
    int m_cached_til_idx = -1;
    #else
    std::vector<std::vector<int>> data;
    #endif
};

#endif
