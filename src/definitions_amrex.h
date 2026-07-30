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

#ifndef DEFINITIONS_AMREX_H_
#define DEFINITIONS_AMREX_H_

#if USE_AMREX
#include <AMReX_MFIter.H>

// =====================================================================
// MFIter_TILING — the single knob for MFIter tiling. Every tile loop in the
// code base (TILE_LOOP and the FIELDLOOP family in looping.h, plus the direct
// MFIter loops in grid_amrex, amrex_solver and printer_CFD) is constructed as
//
//     amrex::MFIter mfi(mf, MFIter_TILING);
//
// so the expression returned below is the only place tiling is decided.
//
// -------------------------------------------------------------------
// Why an object with a conversion, rather than a function or a constant
// -------------------------------------------------------------------
// NOT a function. Written bare — which is how every call site reads — a function
// decays to a function pointer, and that pointer converts to bool -> true. It
// therefore selects MFIter(const FabArrayBase&, bool do_tiling) with do_tiling
// TRUE no matter what the body returns: the knob is inoperative and tiling is
// always on. That is what this replaced. An object does not decay, and the call
// form MFIter_TILING() does not compile, so the knob cannot be bypassed again.
//
// NOT a plain bool constant either, because that would fix the value at static
// initialisation. amrex::TilingIfNotGPU() answers Gpu::notInLaunchRegion(),
// which changes as launch regions open and close, so it has to be re-evaluated
// at each use — hence a conversion operator, which runs per loop.
//
// Overload resolution is unambiguous: reaching either candidate costs the same
// user-defined conversion, after which bool is an identity match and the
// unsigned char of MFIter(const FabArrayBase&, unsigned char) is an integral
// conversion, so the bool overload wins.
// -------------------------------------------------------------------
struct MFIter_Tiling_knob
{
    operator bool() const noexcept
    {
        // ---- change tiling here, and only here ----
        //   true                    : always tile (what the code has been doing)
        //   false                   : one tile per box
        //   amrex::TilingIfNotGPU() : tile on the host, not inside a GPU launch region
        return false;
    }
};

inline constexpr MFIter_Tiling_knob MFIter_TILING{};

// Note on `true`: tiling uses AMReX's default FabArrayBase::mfiter_tile_size =
// (1024000, 8, 8), i.e. boxes split 8 cells deep in y AND z, full extent in x —
// a 32^3 box is 16 tiles. `false` makes amr_tile_lo the box origin and collapses
// the tile-context table to one entry per box, which is a real behaviour change,
// so flip it deliberately rather than in passing.
// =====================================================================

#endif
#endif
