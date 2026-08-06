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
//   2. It can be __restrict-qualified (see view3d<T,true> below). That is the
//      claim that actually unblocks the vector unit, and it is one you can only
//      make per loop, never on a class member.
//
// Both apply to the AMReX build too, where the accessors additionally carry a
// per-access tile-cache branch — see the USE_AMREX section for how the same two
// type names are spelled there.
//
// -------------------------------------------------------------------
// Same names, same call syntax, DIFFERENT layout per build
// -------------------------------------------------------------------
// A converted loop body reads `u(i,j,k)` in both builds. The struct behind it
// is not shared, and deliberately so: this build stores k contiguously
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
    #define REEF3D_VIEW_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
    #define REEF3D_RESTRICT __restrict
    #define REEF3D_VIEW_INLINE __forceinline
#else
    #define REEF3D_RESTRICT
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

#endif
