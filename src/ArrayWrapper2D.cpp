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

#include "ArrayWrapper2D.h"
#include "ArrayWrapper2D_imp.h"
#include "lexer.h"
#include <algorithm>

ArrayWrapper2D::ArrayWrapper2D(lexer *pp) : p(pp)
{
    #if USE_AMREX
    p->register_slice(this);
    #endif
}

ArrayWrapper2D::~ArrayWrapper2D()
{
    #if USE_AMREX
    if(p) p->deregister_slice(this);
    #endif
}

void ArrayWrapper2D::resize(int default_value)
{
    #if USE_AMREX
    m_default = default_value;
    nlevs   = p->nlevs;
    m_ghost = amrex::IntVect(AMREX_D_DECL(p->margin, p->margin, 0));

    m_view.resize(nlevs);
    m_unique.resize(nlevs);

    for(int lev=0; lev<nlevs; ++lev)
        define_level(lev, /*seed_from_coarse=*/false);
    #else
    data.resize(static_cast<std::size_t>(p->imax)*p->jmax, default_value);
    cache_addressing();
    #endif
}

void ArrayWrapper2D::setVal(int val, bool includeGhost)
{
    #if USE_AMREX
    const amrex::IntVect ng = includeGhost ? m_ghost : amrex::IntVect(0);
    for(int lev=0; lev<nlevs; ++lev)
    {
        m_view[lev].setVal(val, 0, 1, ng);
        m_unique[lev].setVal(val, 0, 1, ng);
    }
    #else
    if(includeGhost)
    {
        std::fill(data.begin(), data.end(), val);
    }
    else
    {
        int i,j;
        SLICELOOP4
        {
            operator()(i,j) = val;
        }
    }
    #endif
}

ArrayWrapper2D::operator int *()
{
    #if USE_AMREX
    return m_view[p->level].atLocalIdx(p->amr_local_fab_idx).dataPtr(0);
    #else
    return data.data(); // unshifted: callers index this with IJ, not with i/j
    #endif
}

#if USE_AMREX

// ---------------------------------------------------------------
// fillBoundary — same-level x/y halo exchange.
// Routed through m_unique: reduce the overlapping view onto the
// removeOverlap partition, exchange there, broadcast back. This is what
// makeUnique() already does, so the whole thing is one call.
// ---------------------------------------------------------------
void ArrayWrapper2D::fillBoundary()
{
    for(int lev=0; lev<nlevs; ++lev) makeUnique(lev);
}

// ---------------------------------------------------------------
// inject_from_coarse — piecewise-constant coarse->fine injection.
// No Interpolater (int payload), z ratio pinned to 1 because every level's slab
// sits at z-index 0. Fills fine_mf's valid cells AND the x/y ghost ring inside
// the fine domain (the coarse-fine ring FillBoundary cannot reach); cells
// outside the domain are left alone so their wall flags survive.
// ---------------------------------------------------------------
void ArrayWrapper2D::inject_from_coarse(amrex::iMultiFab& fine_mf,
                                                const amrex::iMultiFab& crse_unique,
                                                int fine_lev)
{
    const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));
    const int ratio_x = ratio2d[0];
    const int ratio_y = ratio2d[1];

    // Coarse data on the fine layout, with enough ghost that the coarsened fine
    // ghost ring is covered. m_ghost is anisotropic (z=0), so drive the widths
    // off m_ghost directly rather than iMultiFab::nGrow() (which assumes
    // isotropic ghosts).
    const int ng_f = m_ghost[0];
    const int rmax = std::max(ratio_x, ratio_y);
    const int ng_c = (ng_f + rmax - 1) / rmax;

    amrex::BoxArray crsnd_fine_ba = amrex::coarsen(fine_mf.boxArray(), ratio2d);
    amrex::iMultiFab crse_on_fine(crsnd_fine_ba, fine_mf.DistributionMap(), 1, ng_c);
    crse_on_fine.ParallelCopy(crse_unique, 0, 0, 1,
                              amrex::IntVect(0), amrex::IntVect(ng_c),
                              p->amrex_geometry[fine_lev-1].periodicity());

    // Domain flattened to the slab: only x/y bounds are meaningful.
    const amrex::Box fine_domain = makeSlab(p->amrex_geometry[fine_lev].Domain(), 2, 0);

    for(amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box fill_box = amrex::grow(mfi.validbox(), m_ghost);
        auto const& fine_arr = fine_mf.array(mfi);
        auto const& crse_arr = crse_on_fine.const_array(mfi);

        amrex::ParallelFor(fill_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(!fine_domain.contains(amrex::IntVect(i,j,k))) return;
            fine_arr(i,j,k,0) = crse_arr(amrex::coarsen(i, ratio_x),
                                         amrex::coarsen(j, ratio_y), k, 0);
        });
    }
}

// ---------------------------------------------------------------
// fillHigherLevels — 2D analogue of ArrayWrapper3D::fillHigherLevels.
// ---------------------------------------------------------------
void ArrayWrapper2D::fillHigherLevels()
{
    for(int lev=1; lev<nlevs; ++lev)
    {
        // Reduce lev-1 first so the source is the authoritative partition.
        makeUnique(lev-1);
        inject_from_coarse(m_view[lev], m_unique[lev-1], lev);
        // Re-reduce so m_unique[lev] agrees with the freshly injected view.
        makeUnique(lev);
    }
}

// ---------------------------------------------------------------
// regrid — rebuild the slab layout, preserving m_unique as the carrier.
// Same MakeNewLevelFromCoarse idiom as slice_amrex::regrid(): seed a new
// level by injection from lev-1, then let retained fine data win.
// ---------------------------------------------------------------
void ArrayWrapper2D::regrid()
{
    nlevs = p->nlevs;
    m_view.resize(nlevs);
    m_unique.resize(nlevs);

    for(int lev=0; lev<nlevs; ++lev)
        define_level(lev, /*seed_from_coarse=*/lev > 0);

    m_cached_level = m_cached_mfi_idx = m_cached_til_idx = -1;
}

void ArrayWrapper2D::define_level(int lev, bool seed_from_coarse)
{
    const auto& dm3d = p->amrex_distribution_mapping[lev];

    // The view layout is the grid's, not ours — sharing it is what makes an
    // MFIter over m_view a valid index into p->slice_owner(lev).
    const amrex::BoxArray& ba2d_view = p->slice_view_boxarray(lev);

    if(m_view[lev].ok()
       && m_view[lev].boxArray() == ba2d_view
       && m_view[lev].DistributionMap() == dm3d)
        return;   // layout unchanged — nothing to rebuild

    amrex::BoxArray ba2d_unique = ba2d_view;
    ba2d_unique.removeOverlap();
    amrex::DistributionMapping dm2d(ba2d_unique);

    amrex::iMultiFab new_unique(ba2d_unique, dm2d, 1, m_ghost);
    new_unique.setVal(m_default);

    // Seed the whole level from the coarser one (MakeNewLevelFromCoarse idiom):
    // a brand-new level, or a newly refined region of an existing one, has no
    // fine data to retain, so without the seed those columns would enter
    // makeUnique() as m_default. Injection from m_unique[lev-1] — already rebuilt
    // for the new layout since define_level runs ascending — stands in for
    // slice_amrex's InterpFromCoarseLevel, which is Real-only.
    if(seed_from_coarse && lev > 0)
        inject_from_coarse(new_unique, m_unique[lev-1], lev);

    new_unique.ParallelCopy(m_unique[lev], 0, 0, 1, 0, 0);   // retained fine data wins
    m_unique[lev] = std::move(new_unique);

    m_view[lev].define(ba2d_view, dm3d, 1, m_ghost);
    m_view[lev].setVal(m_default);

    m_view[lev].ParallelCopy(m_unique[lev], 0, 0, 1, amrex::IntVect(0), m_ghost);
    makeUnique(lev);
}

void ArrayWrapper2D::makeUnique(int lev)
{
    const amrex::iMultiFab& owner = p->slice_owner(lev);

    // Zero non-owner duplicates so the ParallelAdd reduce counts each
    // column exactly once (see slice_amrex::makeUnique).
    for(amrex::MFIter mfi(m_view[lev]); mfi.isValid(); ++mfi)
    {
        auto view = m_view[lev].array(mfi);
        auto own  = owner.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(),
        [=] AMREX_GPU_DEVICE (int i, int j, int) { if(own(i,j,0) == 0) view(i,j,0) = 0; });
    }

    m_unique[lev].setVal(0);
    m_unique[lev].ParallelAdd(m_view[lev]);

    fillHoles(lev);

    m_view[lev].ParallelCopy(m_unique[lev], 0, 0, 1, amrex::IntVect(0), m_ghost,
                             p->amrex_geometry[lev].periodicity());
}

// Columns with no owner at this level (fine level not reaching z-lo) take
// the coarse column's flag by injection — never interpolate a flag.
void ArrayWrapper2D::fillHoles(int lev)
{
    if(lev == 0) return;   // level 0 spans the full column: no holes

    amrex::iMultiFab count(m_unique[lev].boxArray(), m_unique[lev].DistributionMap(), 1, 0);
    count.setVal(0);
    count.ParallelAdd(p->slice_owner(lev));

    const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));
    amrex::BoxArray crsnd_ba = amrex::coarsen(m_unique[lev].boxArray(), ratio2d);
    amrex::iMultiFab crse_on_fine(crsnd_ba, m_unique[lev].DistributionMap(), 1, 0);
    crse_on_fine.ParallelCopy(m_unique[lev-1], 0, 0, 1);

    const int ratio_x = ratio2d[0];
    const int ratio_y = ratio2d[1];

    for(amrex::MFIter mfi(m_unique[lev]); mfi.isValid(); ++mfi)
    {
        auto uni  = m_unique[lev].array(mfi);
        auto cfil = crse_on_fine.const_array(mfi);
        auto num  = count.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if(num(i,j,0) == 0)
                uni(i,j,0) = cfil(amrex::coarsen(i, ratio_x),
                                  amrex::coarsen(j, ratio_y), k, 0);
        });
    }
}

#else

/*!
    * @brief Precomputes the flat addressing used by operator() and operator[].
    *
    * 2D counterpart of ArrayWrapper3D::cache_addressing — see there for the
    * rationale. Here
    *
    *   data[(ii-imin)*jmax + (jj-jmin)]  ==  m_base[ii*m_js + jj]
    *
    * so operator() needs no lexer dereference. imin/jmin are -margin, so m_base
    * lands inside the allocation rather than before it.
*/
void ArrayWrapper2D::cache_addressing() noexcept
{
    m_js   = p->jmax;
    m_base = data.data() - p->imin*m_js - p->jmin;
}
#endif
