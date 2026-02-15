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
    #include <algorithm>
#endif

/*!
 * @brief Base class for field representations in REEF3D.
 *
 * This class provides a generic interface for field operations such as element access,
 * value setting, and boundary condition handling. It is templated to allow for different
 * data types (e.g., double, int).
 */

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

    /*!
     * @brief Accesses an element in the field.
     *
     * Provides a reference to the element at the specified 3D index (ii, jj, kk).
     *
     * @param ii The index in the x-direction.
     * @param jj The index in the y-direction.
     * @param kk The index in the z-direction.
     * @param addOrigin Flag to indicate if the local origin should be added to the indices. Defaults to true.
     * @return T& Reference to the element at the specified location.
     */
#if USE_AMREX
    virtual T& operator()(int ii, int jj, int kk, bool addOrigin=true) noexcept = 0;

    /*!
     * @brief Sets all elements in the field to a specific value.
     *
     * @param val The value to set the field elements to.
     * @param includeGhost Flag to indicate if ghost cells/boundary layers should also be set to this value. Defaults to false.
     */
    virtual void setVal(T val, bool includeGhost = false) = 0;
#else
    inline T& operator()(int ii, int jj, int kk) noexcept {return V[(ii-imin)*jkmax + (jj-jmin)*kmax + kk-kmin];};
    void setVal(T val, bool includeGhost = false) {int i,j,k; if(includeGhost){std::fill(V.begin(),V.end(),val);} else{LOOP{operator()(i,j,k) = val;}}};
#endif

#if USE_AMREX
    /*!
     * @brief Updates the boundary conditions for the field.
     *
     * This method handles the synchronization of ghost cells based valid internal values.
     * Will be deprecated in future versions in favor of AMReX specific boundary exchange methods.
     */
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
