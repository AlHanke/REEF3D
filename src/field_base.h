/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef FIELD_BASE_H_
#define FIELD_BASE_H_

/*!
 * @brief Base class for field representations in REEF3D.
 *
 * This class provides a generic interface for field operations such as element access,
 * value setting, and boundary condition handling. It is templated to allow for different
 * data types (e.g., double, int).
 */

template <typename T>
class field_base
{
public:
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
    virtual T& operator()(int ii, int jj, int kk, bool addOrigin=true) = 0;

    /*!
     * @brief Sets all elements in the field to a specific value.
     *
     * @param val The value to set the field elements to.
     * @param includeGhost Flag to indicate if ghost cells/boundary layers should also be set to this value. Defaults to false.
     */
    virtual void setVal(T val, bool includeGhost = false) = 0;

    /*!
     * @brief Updates the boundary conditions for the field.
     *
     * This method handles the synchronization of ghost cells based valid internal values.
     * Will be deprecated in future versions in favor of AMReX specific boundary exchange methods.
     */
    virtual void FillBoundary() = 0;
};

#endif