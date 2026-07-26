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

#include"wave_lib_hdc.h"

void wave_lib_hdc::allocate_fnpf()
{
    E1.resize(Nx, std::vector<double>(Ny));
    E2.resize(Nx, std::vector<double>(Ny));
    E.resize(Nx, std::vector<double>(Ny));

    F1.resize(Nx, std::vector<double>(Ny));
    F2.resize(Nx, std::vector<double>(Ny));
    F.resize(Nx, std::vector<double>(Ny));
}

void wave_lib_hdc::allocate_cfd()
{
    U1.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    U2.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    U.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    V1.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    V2.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    V.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    W1.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    W2.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));
    W.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    Z.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz)));

    E1.resize(Nx, std::vector<double>(Ny));
    E2.resize(Nx, std::vector<double>(Ny));
    E.resize(Nx, std::vector<double>(Ny));
}
