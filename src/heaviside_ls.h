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

#ifndef HEAVISIDE_LS_H_
#define HEAVISIDE_LS_H_

#include <cmath>

// Heaviside function for level set methods
// Returns 1.0 if phi > eps, 0.0 if phi < -eps, and a smooth transition in between
static inline double heaviside_ls(double phi, double eps)
{
    if (phi > eps)
        return 1.0;
    else if (phi < -eps)
        return 0.0;
    else
        return 0.5 * (1.0 + phi / eps + std::sin(M_PI * phi / eps) / M_PI);
}

#endif