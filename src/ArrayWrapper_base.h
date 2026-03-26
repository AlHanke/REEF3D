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

#ifndef ARRAYWRAPPER_BASE_H_
#define ARRAYWRAPPER_BASE_H_

template <typename T>
class ArrayWrapper_base
{
public:
    virtual ~ArrayWrapper_base() = default;

    virtual void resize(T default_value = 0) = 0;

    virtual T& operator[] (int index) = 0;
    virtual operator T* () = 0;
    virtual void setVal(T val, bool includeGhost = false);

    virtual void fillBoundary() {}; // Only needed for AMReX MultiFabs, can be a no-op for other implementations
    virtual void fillHigherLevels() {}; // Only needed for AMReX MultiFabs, can be a no-op for other implementations
};

#endif
