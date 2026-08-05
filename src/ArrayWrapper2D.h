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

#include <vector>

class lexer;

class ArrayWrapper2D final
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

private:
    /// Recomputes the cached addressing below. Must be called after anything
    /// that resizes @p data. See ArrayWrapper3D::cache_addressing for why the
    /// origin is folded into the pointer and why the stride is a long.
    void cache_addressing() noexcept;

    lexer *p = nullptr;

    std::vector<int> data;

    int *m_base = nullptr; ///< origin-folded base: m_base[ii*m_js + jj]
    long m_js   = 0;       ///< i-stride (jmax); long on purpose
};

#endif
