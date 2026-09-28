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
#include <algorithm>

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

// Antiderivative of heaviside_ls on the band [-eps,eps], G(-eps)=0, G(eps)=eps.
static inline double heaviside_ls_int_band(double s, double eps)
{
    return 0.5 * ((s + eps) + (s*s - eps*eps) / (2.0*eps)
                  - (eps / (M_PI*M_PI)) * (1.0 + std::cos(M_PI * s / eps)));
}

// Mean of heaviside_ls over phi in [a,b], phi linear along the segment. A face density built from
// this makes a column sum of rho_face*dz the exact integral of rho, independent of the grid, so
// every AMR level samples the same continuous hydrostatic profile (the midpoint value
// heaviside_ls(0.5*(a+b)) does not: its error is grid-dependent and shows up as a coarse/fine
// pressure offset above the band). Short intervals (horizontal faces, where a ~ b) use 3-point
// Gauss instead of the difference quotient, which would cancel catastrophically.
static inline double heaviside_ls_avg(double a, double b, double eps)
{
    if (a > b) { const double t = a; a = b; b = t; }
    const double h = b - a;

    if (h < 0.01 * eps)
    {
        const double m = 0.5 * (a + b), r = 0.5 * h * std::sqrt(0.6);
        return (5.0*heaviside_ls(m - r, eps) + 8.0*heaviside_ls(m, eps) + 5.0*heaviside_ls(m + r, eps)) / 18.0;
    }

    double I = 0.0;
    const double lo = std::max(a, -eps), hi = std::min(b, eps);
    if (hi > lo) I += heaviside_ls_int_band(hi, eps) - heaviside_ls_int_band(lo, eps);
    if (b > eps) I += b - std::max(a, eps);
    return I / h;
}

#endif