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
#include <AMReX_Array4.H>

// =====================================================================
// collapse_y_stride — pseudo-2D (knoy==1, j_dir==0) y-index collapse.
//
// The ~2400 legacy stencil sites in src/ index fields as f(i,j+-1..margin,k),
// and every one of them reaches the data through field_amrex / fieldint_amrex
// operator().  In a pseudo-2D run there is only one valid y plane, and the
// physical-BC fillers deliberately skip every y ghost cell (amrex_bc_func.h
// and amrex_bc_func2D.h: "if(!y_dimension_exists && iv[1]!=0)"), so those
// offsets currently resolve to halo memory that nothing ever fills.
//
// Zeroing the Array4's y stride makes every j resolve to the single valid
// plane.  "j*stride.a[0]" is already a runtime multiply in AMReX's addressing
// (AMReX_Array4.H, ptr(i,j,k)), so this costs NOTHING per access — no extra
// instruction, no branch, no clamp.  The whole adjustment lives in the
// cache-refresh cold path, which runs once per FAB/tile change, so a genuine
// 3D run pays one well-predicted not-taken branch per tile and its accessor
// code generation is byte-for-byte unchanged.
//
// Idempotent: once stride.a[0] is 0 the pointer re-anchor is a no-op, so this
// is safe to call from every cache-refresh branch.
// =====================================================================
// =====================================================================
// PSEUDO2D_YGHOST / field_ghost — pseudo-2D y-ghost reduction.
//
// With knoy==1 and margin==3, an isotropic ghost width gives every field a y
// extent of 1+2*3 = 7 planes of which ONE is valid: ~86% of each fab is halo
// that the physical-BC fillers deliberately never write.  Narrowing the y ghost
// shrinks the fab, the FillBoundary copy volume and the FillPatch interpolation
// work by the same factor.
//
// Reading past the valid plane stays safe regardless of this width because
// collapse_y_stride() below folds every j onto it.  The value here is therefore
// governed by what AMReX's OWN y-stencils need, not by REEF3D's:
//
//   1 (default) — 7 planes -> 3, a 2.33x reduction.  Keeps the single coarse
//                 ghost that cell_cons_interp's slope stencil reads inside
//                 FillPatchTwoLevels / InterpFromCoarseLevel (grid_amrex.cpp
//                 :200, :289, :1121 and field_amrex.h :750, :897), so the
//                 coarse-fine path behaves exactly as it does today.
//   0           — 7 planes -> 1, the full 7x.  Unconditionally safe when
//                 nlevs==1 (no interpolation happens at all).  For nlevs>1 it
//                 must be validated against the C-F interp path first: the
//                 cell_cons_interp slope stencil would then read a coarse y
//                 ghost that no longer exists.
//
// A 3D run is unaffected: j_dir==1 returns the isotropic margin unchanged.
// =====================================================================
inline constexpr int PSEUDO2D_YGHOST = 1;

/// Ghost width for 3D field storage — full margin in x/z, narrowed in y for a
/// pseudo-2D run. Every field_amrex / fieldint_amrex allocation goes through this.
AMREX_FORCE_INLINE amrex::IntVect field_ghost(int margin, bool j_dir) noexcept
{
    return amrex::IntVect(margin, j_dir ? margin : PSEUDO2D_YGHOST, margin);
}

/// Ghost width for the horizontal slice family (slice_amrex, sliceint_amrex,
/// ArrayWrapper2D). These live in 3D fabs with a single z plane, so their y is
/// component 1 exactly as for the volume fields and collapse_y_stride() applies
/// to them unchanged — only the z component differs, which is always 0 here.
AMREX_FORCE_INLINE amrex::IntVect slice_ghost(int margin, bool j_dir) noexcept
{
    return amrex::IntVect(AMREX_D_DECL(margin, j_dir ? margin : PSEUDO2D_YGHOST, 0));
}

template <typename T>
AMREX_FORCE_INLINE void collapse_y_stride(amrex::Array4<T>& arr, int margin) noexcept
{
    // Re-anchor the base pointer onto y == 0 BEFORE dropping the stride: "p" is
    // anchored at "begin", so with a zero y stride every access would otherwise
    // land on the lowest ghost plane (y == -margin) instead of the valid one.
    arr.p -= arr.begin.vect[1] * arr.stride.a[0];
    arr.stride.a[0] = 0;

    // With the y stride at 0 the y term drops out of both addressing paths in
    // AMReX_Array4.H (release "idx1-idx0" and debug "(j-begin)*stride"), so
    // begin/end in y now only feed AMREX_ARRAY4_INDEX_ASSERT.  Widen them to the
    // halo range the legacy stencils actually reach so debug builds stay quiet.
    arr.begin.vect[1] = -margin;
    arr.end.vect[1]   =  margin + 1;
}

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
        return amrex::TilingIfNotGPU();
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
