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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef FLUXLIM_H_
#define FLUXLIM_H_

#include"increment.h"

class lexer;

class fluxlim
{
public:
    template<typename GernericField>
    inline double iphi(const GernericField& b, int n1, int n2, int q1, int q2)
    {
        const int i=increment::i, j=increment::j, k=increment::k;
        return phi_impl(b(i+n1,j,k), b(i+n2,j,k), b(i+q1,j,k), b(i+q2,j,k), b(i,j,k));
    }

    template<typename GernericField>
    inline double jphi(const GernericField& b, int n1, int n2, int q1, int q2)
    {
        const int i=increment::i, j=increment::j, k=increment::k;
        return phi_impl(b(i,j+n1,k), b(i,j+n2,k), b(i,j+q1,k), b(i,j+q2,k), b(i,j,k));
    }

    template<typename GernericField>
    inline double kphi(const GernericField& b, int n1, int n2, int q1, int q2)
    {
        const int i=increment::i, j=increment::j, k=increment::k;
        return phi_impl(b(i,j,k+n1), b(i,j,k+n2), b(i,j,k+q1), b(i,j,k+q2), b(i,j,k));
    }

protected:
    virtual inline double phi_impl(double vn1, double vn2, double vq1, double vq2, double vcell) = 0;
};

#endif
