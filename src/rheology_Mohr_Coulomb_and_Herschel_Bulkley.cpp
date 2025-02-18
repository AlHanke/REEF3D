/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>
#include<math.h>

// https://doi.org/10.5194/gmd-9-2909-2016
// Three phase viscosity - air, slurry and gravel
// Air phase is covered in fluid_update_rheology.cpp
double rheology_f::Mohr_Coulomb_and_Herschel_Bulkley(lexer* p, fdm* a, ghostcell* pgc)
{
    double shear_rate = sqrt(2)*strainterm(p,a);
    const double b = 0.3;
    double n = 0; // calibration : shear thinning vs shear thickening
    const double C = 0.0; // volumetric solid concentration: the volume of all solid particles including fine material relative to the volume of the debris-flow material including water m^3/m^3
    const double tau_00 = 0.0; // calibration : grid resolution sensitivity of the shear rate
    double tau_0 = C<0.47?tau_00:tau_00*exp(5*(C-0.47));
    const double C_kaolinite_chlorite = 0.0;
    const double C_illite = 0.0;
    const double C_montmorillonite = 0.0;
    double P_0 = C_kaolinite_chlorite + 1.3*C_illite + 1.7*C_montmorillonite;
    double P_1 = P_0>0.27 ? 0.7*P_0 : P_0;
    double tau_y = tau_0*C*C*exp(22*C*P_1);
    double k = b*tau_y; // d_max<0.4mm
    double mu_2 = k * pow(fabs(shear_rate),n-1) + tau_y*pow(fabs(shear_rate),-1);
    const double mu_0 = 0.0;
    mu_2 = std::min(mu_2,mu_0);
    const double rho_2 = 0.0;
    double ny_2 = mu_2/rho_2;

    const double mu_min = 0.0; // minimal dynamic viscosity
    const double delta = 0.0; // internal friction angle
    const double m_y = 0.2; // ]0,1]
    double mu_3 = mu_min + a->press(i,j,k)*sin(delta)/shear_rate*(1-exp(-m_y*shear_rate));
    // ToDo: limiter
    const double rho_3 = 0.0;
    double ny_3 = mu_3/rho_3;

    double a_2 = 0.0;
    double a_3 = 0.0;
    return a_2*ny_2+a_3*ny_3;
}