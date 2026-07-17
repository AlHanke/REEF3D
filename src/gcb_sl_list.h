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

#ifndef GCB_SL_LIST_H_
#define GCB_SL_LIST_H_

// =====================================================================
// Ghost-cell boundary list entries.
//
// A gcb list used to be an int** with a hand-managed row count: gcbsl1[n][2]
// was the face direction, gcbsl1[n][3] the boundary condition, and knowing
// that was the reader's problem. The entry types below give those columns
// names, so the same std::vector serves every dimensionality and usecase:
//
//   gcb_sl_cs_bc_list = std::vector<gcb_sl_cs_bc>  2D slice lists: gcbsl1/2/4
//   gcb_list          = std::vector<gcb_field>     3D field lists: gcb1/2/3/4
//
//   gcb_sl_list     = std::vector<gcb_sl>     2D slice inflow list: gcbslin
//   gcb_sl_cs_list  = std::vector<gcb_sl_cs>  2D slice outflow list: gcbslout
//
// Nothing inspects the entry generically, so the only thing that varies is
// the payload (gcb_field carries k, gcb_sl does not). Loop macros and any
// code that only reads .i/.j/.cs/.bc work against either.
//
// This header must not reference lexer: it is included into lexer.h before
// the lexer class is formed.
// =====================================================================

#include <vector>

// One slice inflow ghost-cell boundary entry. Replaces the old int[2] row,
// whose [0]/[1] were i/j.
struct gcb_sl
{
    int i = 0;
    int j = 0;
};

// One slice outflow ghost-cell boundary entry. Replaces the old int[3] row,
// whose [0]/[1] were i/j/cs.
struct gcb_sl_cs
{
    int i  = 0;
    int j  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4
};

// One slice ghost-cell boundary entry. Replaces the old int[4] row, whose
// [0]/[1]/[2]/[3] were i/j/cs/bc.
struct gcb_sl_cs_bc
{
    int i  = 0;
    int j  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4
    int bc = 21;   ///< boundary condition flag; 21 = default wall
};

// One field ghost-cell boundary entry. Replaces the old int[5] row, whose
// [0]/[1]/[2]/[3]/[4] were i/j/k/cs/bc.
struct gcb_field
{
    int i  = 0;
    int j  = 0;
    int k  = 0;
    int cs = 0;    ///< face direction: X_NEG=1, Y_POS=2, Y_NEG=3, X_POS=4, Z_POS=5, Z_NEG=6
    int bc = 21;   ///< boundary condition flag; 21 = default wall
};

// Number of entries, as int — the loop macros compare against an int
// counter, and the dev target builds with -Wsign-conversion.
template<typename ENTRY>
inline int gcb_ssize(const std::vector<ENTRY> &list) noexcept
{
    return static_cast<int>(list.size());
}

using gcb_sl_cs_bc_list = std::vector<gcb_sl_cs_bc>;   ///< 2D: gcbsl1/2/4
using gcb_list          = std::vector<gcb_field>;      ///< 3D: gcb1/2/3/4

using gcb_sl_list    = std::vector<gcb_sl>;      ///< 2D: gcbslin
using gcb_sl_cs_list = std::vector<gcb_sl_cs>;   ///< 2D: gcbslout

#endif
