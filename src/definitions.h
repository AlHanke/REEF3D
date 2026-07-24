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

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

enum {
    WATER_FLAG = 10,
    AIR_FLAG = -1,
    INFLOW_FLAG = -3,
    OUTFLOW_FLAG = -4,
    FLT_FLAG = -17,
    TOPO_FLAG = -18,
    SOLID_FLAG = -19,
    OBJ_FLAG = -20
};

enum {X_NEG=1, X_POS=4, Y_NEG=3, Y_POS=2, Z_NEG=5, Z_POS=6};

enum {BC_INFLOW=1, BC_OUTFLOW=2, BC_SYMMETRY=3, BC_WAVEGEN=6, BC_NUMBEACH=7, BC_WALL=21};

inline constexpr double PI = 3.14159265359;
inline constexpr double EE = 2.71828182846;

/// Storage layout of a field, and the `location` value carried by
/// lexer::mf_registry / imf_registry entries (field_amrex passes
/// static_cast<int>(DataLocation)).
///
/// CELL_CENTERED and FACE_X/Y/Z all live in a cell-centred data structure — the
/// staggered ones by index convention, re-staggered only at interpolation
/// time.
enum class DataLocation : unsigned int { CELL_CENTERED = 0, FACE_X = 1, FACE_Y = 2, FACE_Z = 3, NODE_Z = 7 };

#endif
