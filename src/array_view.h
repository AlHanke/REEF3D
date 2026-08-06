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

#ifndef ARRAY_VIEW_H_
#define ARRAY_VIEW_H_

#include <type_traits>

// =====================================================================
// view3d / view2d — capture-by-value windows onto a container's storage.
//
// A view is a raw pointer plus the strides needed to reach an element, read
// ONCE and then carried in registers, instead of reached through `this` on
// every access. It is what a loop body should hold: materialise the view where
// the container is known (outside the nest), index it by value inside.
//
// The point is not the addressing arithmetic — field_base and ArrayWrapper3D
// already fold origin and strides into a member pointer, so a single access
// costs the same either way (see Accessor_Performance.md). The point is what a
// view lets you say that a member cannot:
//
//   1. Nothing in the loop body can invalidate it. A view is a local whose
//      address is never taken, so no store and no opaque call can force the
//      base and strides to be reloaded.
//
//   2. In the plain build it can be __restrict-qualified (view3d<T,true>). That
//      is the claim that actually unblocks the vector unit, and it is one you
//      can only make per loop, never on a class member. In the AMReX build the
//      claim is already made for you — amrex::Array4 declares its base pointer
//      `T* AMREX_RESTRICT p` — so there the parameter is inert and only item 1
//      is new.
//
// Both apply harder in the AMReX build, where operator() additionally runs
// refresh_cache_if_needed() on EVERY access: three loads out of lexer, two
// compares, a branch, and conditional stores into mutable members. Those stores
// are a hard vectorisation blocker — the compiler cannot prove the branch is
// never taken, so the body carries loop-carried state and cannot be
// if-converted. Hoisting a view per tile removes all of it. FIELDLOOP already
// does exactly this for the double fields via const_array(_fl_mfi); the int flag
// wrappers behind PCHECK/UCHECK/... are the ones that never got the treatment.
//
// -------------------------------------------------------------------
// Same names, same call syntax, DIFFERENT layout per build
// -------------------------------------------------------------------
// A converted loop body reads `u(i,j,k)` in both builds. The struct behind it
// is not shared, and deliberately so: the plain build stores k contiguously
// (index = i*is + j*js + k), amrex::Array4 stores i contiguously
// (index = i + j*js + k*ks). Forcing one layout on both would cost a runtime
// multiply against a stride that is 1 by construction in whichever build lost.
// So each build folds its own unit stride into the indexing expression and
// keeps exactly two runtime strides.
//
// Consequence: views are NOT aggregate-initialised by callers. Build them with
// the make_* factories below, whose signatures differ between builds; the
// containers that call them are already per-build classes (field_base vs
// field_amrex, the two halves of ArrayWrapper3D), so nothing new is #if'd.
//
// -------------------------------------------------------------------
// Indices are REEF3D-local in both builds
// -------------------------------------------------------------------
// (i,j,k) are what IJKLOOP produces — 0..IMAX_LOOP, ghost cells reached with
// negative values — never AMReX global indices. That is what the containers'
// operator() already takes in both builds (the plain one folds against
// imin/jmin/kmin, the AMReX one adds m_cached_ox/oy/oz), and it is what
// LocalArr4 in field.h takes. The AMReX factories therefore fold the tile
// origin into the base, exactly as LocalArr4 subtracts it per access.
// =====================================================================

// ---------------------------------------------------------------------
// Stride width. Long for every payload, in every container that hands out a
// view — see field_base::stride_t for the measured rationale. The short of it:
// an int store may alias an int stride member under TBAA, and every container a
// loop body touches must agree on the width or their index computations cannot
// be shared. This is the single definition; field_base and slice_base alias it.
// ---------------------------------------------------------------------
using stride_t = long;

#if defined(__GNUC__) || defined(__clang__)
    #define REEF3D_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define REEF3D_RESTRICT __restrict
#else
    #define REEF3D_RESTRICT
#endif

#if USE_AMREX
    // Views are captured by value into ParallelFor lambdas, so the accessor has
    // to be callable on device.
    #include <AMReX_Array4.H>
    #include <AMReX_Dim3.H>
    #define REEF3D_VIEW_INLINE AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
#elif defined(__GNUC__) || defined(__clang__)
    #define REEF3D_VIEW_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
    #define REEF3D_VIEW_INLINE __forceinline
#else
    #define REEF3D_VIEW_INLINE inline
#endif

namespace view_detail
{
    // Lazy conditional: picking between two wrappers and then taking ::type
    // instantiates only the selected one, which keeps the restrict qualifier
    // attached to the typedef rather than being stripped by a direct
    // conditional_t<R, T* REEF3D_RESTRICT, T*>.
    template<typename T> struct plain_ptr    { using type = T*; };
    template<typename T> struct restrict_ptr { using type = T* REEF3D_RESTRICT; };

    template<typename T, bool R>
    using ptr_t = typename std::conditional_t<R, restrict_ptr<T>, plain_ptr<T>>::type;
}

// =====================================================================
// Restrict, and when you may ask for it
// =====================================================================
// view3d<T,true> asserts that, for as long as the view is live, the memory it
// covers is reached ONLY through this view. Violating that is undefined
// behaviour with no diagnostic — the loop silently computes the wrong answer
// under -O3 and the right one under -O0.
//
// The parameter only does anything in the plain build; amrex::Array4 already
// restrict-qualifies its base pointer, so on that path the rules below describe
// a claim AMReX is making on your behalf whether you opt in or not.
//
// Safe: distinct field objects in one body, which is the common shape.
//     FIELDLOOP(frk1, FIELD_CONST(ls); FIELD_CONST_MEMBER(a,L), ...)
//
// NOT safe, and it occurs in this code base:
//   - the same container appearing as both the mutated and a const operand,
//     e.g. FIELDLOOP(a->u, FIELD_CONST_MEMBER(a,u), ...). Two views, one
//     allocation.
//   - staggered neighbours of the same field (u(i,j,k) and u(i+1,j,k)) are
//     fine — that is one view, not two.
//   - anything reached through a second path in the body: a ghostcell call, a
//     legacy helper taking the field by reference, ipol/aij on the same field.
//
// Default is false. Opt in per loop, once you have checked the operand list.
// =====================================================================

#if !USE_AMREX

/*!
 * @brief 3D window: `base[i*is + j*js + k]`, k contiguous.
 *
 * Origin-folded — the caller's (i,j,k) are the same indices operator() takes,
 * ghost cells included, so the fold has already subtracted imin/jmin/kmin.
 * Build with make_view3d() / fold_view3d().
 *
 * @tparam T Element type. Const-ness is carried in T: a read-only view is
 *           view3d<const double>, not a separate ConstView type.
 * @tparam Restrict See the note above. Do not set it without checking the
 *           operand list of the loop it is used in.
 */
template<typename T, bool Restrict = false>
struct view3d
{
    view_detail::ptr_t<T,Restrict> base;
    stride_t is;   ///< i-stride (jmax*kmax)
    stride_t js;   ///< j-stride (kmax); k has unit stride and is folded in

    REEF3D_VIEW_INLINE T& operator()(int i, int j, int k) const noexcept
    { return base[i*is + j*js + k]; }

    /// Read-only rebind. Lets a view3d<T> bind where a view3d<const T> is
    /// wanted without the caller restating the strides. The target type is
    /// spelled against the template parameter U, not T, so that when T is
    /// already const the declaration stays dependent and is never formed —
    /// spelling it view3d<const T,...> would instantiate a conversion to the
    /// enclosing type itself and trip -Wclass-conversion under the dev target.
    template<typename U = T, typename = std::enable_if_t<!std::is_const_v<U>>>
    REEF3D_VIEW_INLINE operator view3d<const U,Restrict>() const noexcept
    { return view3d<const U,Restrict>{base, is, js}; }
};

/*!
 * @brief 2D window: `base[i*is + j]`, j contiguous. 2D analogue of view3d.
 */
template<typename T, bool Restrict = false>
struct view2d
{
    view_detail::ptr_t<T,Restrict> base;
    stride_t is;   ///< i-stride (jmax); j has unit stride and is folded in

    REEF3D_VIEW_INLINE T& operator()(int i, int j) const noexcept
    { return base[i*is + j]; }

    /// Read-only rebind; see view3d's, including why the target names U.
    template<typename U = T, typename = std::enable_if_t<!std::is_const_v<U>>>
    REEF3D_VIEW_INLINE operator view2d<const U,Restrict>() const noexcept
    { return view2d<const U,Restrict>{base, is}; }
};

// ---------------------------------------------------------------------
// Factories
//
// make_* take an ALREADY-FOLDED base — the m_base every container caches in its
// cache_addressing() — and just pair it with the strides.
//
// fold_* do the fold themselves, from the unfolded allocation and the grid
// minima. This is the arithmetic currently copy-pasted into four separate
// cache_addressing() bodies (field_base, slice_base, ArrayWrapper3D,
// ArrayWrapper2D); those can be reduced to a fold_* call so the
// identity `V[(i-imin)*is + (j-jmin)*js + k-kmin] == base[i*is + j*js + k]` is
// stated once.
// ---------------------------------------------------------------------

template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view3d<T,Restrict> make_view3d(T* folded_base, stride_t is, stride_t js) noexcept
{ return view3d<T,Restrict>{folded_base, is, js}; }

template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view2d<T,Restrict> make_view2d(T* folded_base, stride_t is) noexcept
{ return view2d<T,Restrict>{folded_base, is}; }

/// @param data Start of the allocation (index (imin,jmin,kmin)).
template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view3d<T,Restrict> fold_view3d(T* data, int imin, int jmin, int kmin,
                                                  stride_t is, stride_t js) noexcept
{
    return view3d<T,Restrict>{data - static_cast<stride_t>(imin)*is
                                   - static_cast<stride_t>(jmin)*js - kmin,
                              is, js};
}

/// @param data Start of the allocation (index (imin,jmin)).
template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view2d<T,Restrict> fold_view2d(T* data, int imin, int jmin, stride_t is) noexcept
{
    return view2d<T,Restrict>{data - static_cast<stride_t>(imin)*is - jmin, is};
}

#else // USE_AMREX

// =====================================================================
// AMReX build: the same two names over an amrex::Array4 and a tile origin.
//
// This is field.h's LocalArr4 / LocalArr4Const / LocalArr4IntConst, made into
// one template. Those three are the right idea and the right mechanism — they
// are already what FIELDLOOP_INC hoists per tile — but they cover three cells
// of a 2D/3D x const/mutable x double/int matrix that has twelve, and the three
// disagree with each other: LocalArr4 returns double&, LocalArr4Const returns
// double BY VALUE, so a const view cannot stand in for a mutable one and
// const-ness cannot be a property of the element type. Writing the missing nine
// by hand is the proliferation this file exists to avoid; the intended end state
// is that the three in field.h become aliases for view3d and then disappear.
//
// Deliberately NOT folded to a raw base pointer here, unlike the plain build.
// Folding would save one add per index component per access, and would cost the
// AMREX_ARRAY4_INDEX_ASSERT bounds check that Array4::operator() carries under
// AMREX_DEBUG. That is a bad trade: the win that matters on this path is
// hoisting the access out of ArrayWrapper3D::refresh_cache_if_needed at all,
// which delegating to Array4 already achieves in full, and the offsets hoist out
// of the nest anyway once the view is a by-value local.
// =====================================================================

/*!
 * @brief 3D window over one tile of an amrex::Array4.
 *
 * Indices are REEF3D-local (see the header note): `v(i,j,k)` reaches
 * `arr(i+ox, j+oy, k+oz, comp)`, so the same body text works in both builds.
 *
 * @tparam T Element type; const-ness rides here, as in the plain build.
 * @tparam Restrict Accepted so that call sites are identical across builds, but
 *         inert: amrex::Array4 already declares its base pointer
 *         `T* AMREX_RESTRICT p`, so the claim is made for us and there is
 *         nothing left to opt into. The aliasing caveats in the note above still
 *         describe what AMReX is asserting on your behalf.
 */
template<typename T, bool Restrict = false>
struct view3d
{
    amrex::Array4<T> arr;
    int ox, oy, oz;   ///< tile origin = p->amr_tile_lo.{x,y,z}
    int comp;         ///< component, for the shared-MultiFab views

    REEF3D_VIEW_INLINE T& operator()(int i, int j, int k) const noexcept
    { return arr(i + ox, j + oy, k + oz, comp); }

    /// Read-only rebind; see the plain build's, including why the target names
    /// U. Array4<T> -> Array4<const T> is AMReX's own implicit conversion.
    template<typename U = T, typename = std::enable_if_t<!std::is_const_v<U>>>
    REEF3D_VIEW_INLINE operator view3d<const U,Restrict>() const noexcept
    { return view3d<const U,Restrict>{arr, ox, oy, oz, comp}; }
};

/*!
 * @brief 2D window over one tile of a z-slab Array4 (slice_amrex,
 *        ArrayWrapper2D). The slab's single k index is pinned to 0,
 *        matching slice_amrex::operator().
 */
template<typename T, bool Restrict = false>
struct view2d
{
    amrex::Array4<T> arr;
    int ox, oy;
    int comp;

    REEF3D_VIEW_INLINE T& operator()(int i, int j) const noexcept
    { return arr(i + ox, j + oy, 0, comp); }

    template<typename U = T, typename = std::enable_if_t<!std::is_const_v<U>>>
    REEF3D_VIEW_INLINE operator view2d<const U,Restrict>() const noexcept
    { return view2d<const U,Restrict>{arr, ox, oy, comp}; }
};

// ---------------------------------------------------------------------
// Factories. The plain build's fold_* have no counterpart here — there is no
// unfolded allocation to fold, the tile origin is the whole of the offset.
//
// Callers pass p->amr_tile_lo.{x,y,z}, which is exactly what the containers'
// refresh_cache_if_needed() caches into m_cached_ox/oy/oz and what
// _FIELDLOOP_INC_IMPL computes as ox/oy/oz from the tilebox smallEnd.
// ---------------------------------------------------------------------

template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view3d<T,Restrict> make_view3d(amrex::Array4<T> const& arr,
                                                  int ox, int oy, int oz, int comp = 0) noexcept
{ return view3d<T,Restrict>{arr, ox, oy, oz, comp}; }

template<bool Restrict = false, typename T>
REEF3D_VIEW_INLINE view2d<T,Restrict> make_view2d(amrex::Array4<T> const& arr,
                                                  int ox, int oy, int comp = 0) noexcept
{ return view2d<T,Restrict>{arr, ox, oy, comp}; }

#endif // USE_AMREX

#endif
