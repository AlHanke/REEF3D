/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef SLICEINT_AMREX_H_
#define SLICEINT_AMREX_H_

#include "sliceint.h"
#include "slice_amrex_prim.h"
#include "lexer.h"
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Reduce.H>
#include <AMReX_Vector.H>

/// Integer counterpart of slice_amrex: an amrex::iMultiFab-backed 2D slab
/// carrying the same view / unique / owner / ksurf machinery.
///
/// Layout is identical to slice_amrex's — each 3D box makeSlab-flattened to
/// z-index 0 on the 3D DistributionMapping — so the MFIter LocalIndex that
/// TILE_LOOP publishes addresses this container too.
///
/// Two deliberate divergences from slice_amrex, both forced by the payload
/// being int rather than amrex::Real:
///
///   * m_ksurf is an iMultiFab here, not a MultiFab. It holds a k cell index,
///     which is an integer, and buildOwnerMask only ever compares it against
///     the integer z-bounds of the un-flattened box.
///   * the coarse fill is hand-rolled (interp_from_coarse) instead of going
///     through amrex::InterpFromCoarseLevel with &amrex::pc_interp. AMReX has
///     no interpolater that accepts integer data: Interpolater::interp takes
///     FArrayBox and MFInterpolater::interp takes MultiFab, and
///     PhysBCFunctNoOp's operator() is MultiFab-only as well. See
///     interp_from_coarse for why injection reproduces pc_interp exactly.
///
/// There is no BCRec / FillDomainBoundary path: that machinery exists on
/// slice_amrex to evaluate physical boundary conditions on a floating-point
/// field, and an integer slice carries flags rather than a solved quantity.
class sliceint_amrex : public slice_amrex_prim, public sliceint
{
public:
    sliceint_amrex(lexer *pp)
    {
        p = pp;
        p->register_slice(this);
        nlevs = p->nlevs;
        m_ghost = amrex::IntVect(AMREX_D_DECL(p->margin, p->margin, 0));

        m_view.resize(nlevs);
        m_unique.resize(nlevs);
        m_ksurf.resize(nlevs);
        m_owner.resize(nlevs);
        m_ba = &p->amrex_box_array;

        for (int lev = 0; lev < nlevs; ++lev)
        {
            const auto& ba3d = (*m_ba)[lev];
            const auto& dm3d = p->amrex_distribution_mapping[lev];

            amrex::BoxList bl;
            for(int i=0; i<ba3d.size(); ++i)
            {
                bl.push_back(makeSlab(ba3d[i], 2, 0));
            }

            amrex::BoxArray ba2d_view(std::move(bl));
            m_view[lev].define(ba2d_view, dm3d, 1, m_ghost);
            m_view[lev].setVal(0, 0, 1, m_ghost);
            m_owner[lev].define(ba2d_view, dm3d, 1, 0);
            m_owner[lev].setVal(0, 0, 1, 0);

            amrex::BoxArray ba2d_unique = ba2d_view;
            ba2d_unique.removeOverlap();
            amrex::DistributionMapping dm2d(ba2d_unique);
            m_unique[lev].define(ba2d_unique, dm2d, 1, m_ghost);
            m_unique[lev].setVal(0, 0, 1, m_ghost);
            // ksurf lives in the VIEW layout: buildOwnerMask reads it under an
            // MFIter over m_view, and each overlapping view box needs its own copy.
            m_ksurf[lev].define(ba2d_view, dm3d, 1, m_ghost);
            m_ksurf[lev].setVal(0, 0, 1, m_ghost);

            buildOwnerMask(lev);
            makeUnique(lev);
            checkInvariants(lev);
        }
    };
    virtual ~sliceint_amrex()
    {
        if(p)
        {
            p->deregister_slice(this);
        }
    };

    sliceint_amrex(const sliceint_amrex&) = delete;
    sliceint_amrex& operator=(const sliceint_amrex&) = delete;
    sliceint_amrex(sliceint_amrex&&) = delete;
    sliceint_amrex& operator=(sliceint_amrex&&) = delete;

    inline int& operator()(int ii, int jj) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, 0, 0);
    };

    inline const int& operator()(int ii, int jj) const noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, 0, 0);
    };

    void regrid() override final
    {
        nlevs = p->nlevs;
        m_view.resize(nlevs);
        m_unique.resize(nlevs);
        m_ksurf.resize(nlevs);
        m_owner.resize(nlevs);

        for (int lev = 0; lev < nlevs; ++lev)
        {
            const auto& ba3d = (*m_ba)[lev];
            const auto& dm3d = p->amrex_distribution_mapping[lev];

            amrex::BoxList bl;
            for(int i=0; i<ba3d.size(); ++i)
            {
                bl.push_back(makeSlab(ba3d[i], 2, 0));
            }

            amrex::BoxArray ba2d_view(std::move(bl));
            if(m_view[lev].boxArray() != ba2d_view || m_view[lev].DistributionMap() != dm3d)
            {
                amrex::BoxArray ba2d_unique = ba2d_view;
                ba2d_unique.removeOverlap();
                amrex::DistributionMapping dm2d(ba2d_unique);

                // Rebuild the authoritative unique data for the new layout, following
                // the AmrCore RemakeLevel/MakeNewLevelFromCoarse idiom: seed the whole
                // level from the coarser one, then let retained fine data win wherever
                // the previous layout had it.
                //
                // The coarse seed is not optional: a brand-new level (and any newly
                // refined region of an existing level) has no fine data to preserve, so
                // without it those columns would enter makeUnique() as zeros. fillHoles()
                // would not rescue them either — it only fills columns with NO owner,
                // and the owned columns are exactly the ones straddling the surface.
                amrex::iMultiFab new_unique(ba2d_unique, dm2d, 1, m_ghost);
                new_unique.setVal(0, 0, 1, m_ghost);

                if (lev > 0)
                {
                    interp_from_coarse(new_unique, lev);
                }

                // Retained fine data overwrites the coarse seed where it exists.
                // (No-op for a brand-new level: ParallelCopy early-outs on an empty src.)
                new_unique.ParallelCopy(m_unique[lev], 0, 0, 1, 0, 0);
                m_unique[lev] = std::move(new_unique);

                // View, owner and ksurf are caches rebuilt from scratch for the new layout.
                m_view[lev].define(ba2d_view, dm3d, 1, m_ghost);
                m_view[lev].setVal(0, 0, 1, m_ghost);
                m_owner[lev].define(ba2d_view, dm3d, 1, 0);
                m_owner[lev].setVal(0, 0, 1, 0);
                m_ksurf[lev].define(ba2d_view, dm3d, 1, m_ghost);
                m_ksurf[lev].setVal(0, 0, 1, m_ghost);

                // Seed the fresh view from the preserved unique data, then rebuild
                // the mask and re-reduce so view/unique are consistent again.
                m_view[lev].ParallelCopy(m_unique[lev], 0, 0, 1, amrex::IntVect(0), m_ghost);
                buildOwnerMask(lev);
                makeUnique(lev);
                checkInvariants(lev);
            }
        }

        m_cached_level = -1;
        m_cached_mfi_idx = -1;
        m_cached_til_idx = -1;
    };

private:
    void buildOwnerMask(int lev)
    {
        const auto& ba3d = (*m_ba)[lev];
        for (amrex::MFIter mfi(m_view[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& b3 = ba3d[mfi.index()];               // the un-flattened box
            auto own  = m_owner[lev].array(mfi);
            auto ksurf = m_ksurf[lev].const_array(mfi);
            // ksurf(i,j): surface cell index from phi/eta (injected from coarse right after regrid)
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int) {
                own(i,j,0) = (b3.smallEnd(2) <= ksurf(i,j,0) && ksurf(i,j,0) <= b3.bigEnd(2));
            });
        }
    }

    void makeUnique(int lev)
    {
        // Zero every non-owner duplicate first: the broadcast at the end of the
        // previous call filled them with the column value, so without this the
        // ParallelAdd below would sum each column more than once. After masking,
        // exactly one view cell per column is nonzero and the reduce is exact.
        for (amrex::MFIter mfi(m_view[lev]); mfi.isValid(); ++mfi)
        {
            auto view = m_view[lev].array(mfi);
            auto own  = m_owner[lev].const_array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int) {
                if(own(i,j,0) == 0) view(i,j,0) = 0;
            });
        }

        m_unique[lev].setVal(0);
        m_unique[lev].ParallelAdd(m_view[lev]);

        // Fill columns with no owner at this level (holes) from the coarser level.
        fillHoles(lev);

        m_view[lev].ParallelCopy(m_unique[lev], 0, 0, 1, amrex::IntVect(0), m_ghost, p->amrex_geometry[lev].periodicity());
    }

    void fillHoles(int lev)
    {
        if(lev == 0) return;   // level 0 spans the full column: no coarser source, no holes

        // Hole mask in the unique layout: columns with no owner at this level.
        amrex::iMultiFab count(m_unique[lev].boxArray(), m_unique[lev].DistributionMap(), 1, 0);
        count.setVal(0);
        count.ParallelAdd(m_owner[lev]);

        // Inject the coarser slice onto this level's unique layout.
        amrex::iMultiFab coarse_fill(m_unique[lev].boxArray(),
                                     m_unique[lev].DistributionMap(), 1, 0);
        coarse_fill.setVal(0);
        interp_from_coarse(coarse_fill, lev);

        // Blend: keep reduced (owned) values; take the injected coarse value in holes.
        for (amrex::MFIter mfi(m_unique[lev]); mfi.isValid(); ++mfi)
        {
            auto uni  = m_unique[lev].array(mfi);
            auto cfil = coarse_fill.const_array(mfi);
            auto num  = count.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int) {
                if(num(i,j,0) == 0) uni(i,j,0) = cfil(i,j,0);
            });
        }
    }

    /// Piecewise-constant (injection) fill of @p dst — valid cells plus its own
    /// ghosts — from m_unique[lev-1].
    ///
    /// slice_amrex reaches for amrex::InterpFromCoarseLevel with &amrex::pc_interp
    /// here; that is unavailable for integer data, because Interpolater::interp
    /// takes an FArrayBox, MFInterpolater::interp takes a MultiFab, and the
    /// PhysBCFunctNoOp those calls need is MultiFab-only too. Nothing is lost:
    /// pc_interp's stencil is a single coarse cell, so it IS this loop.
    ///
    /// Doing it by hand also removes the hazard slice_amrex has to warn about.
    /// A slope-based interpolater grows its coarse patch by one cell in every
    /// direction and would read coarse cells at z = +-1 that the slab does not
    /// have (m_ghost[2] == 0) — uninitialised reads that a z-ratio of 1 cannot
    /// neutralise. This loop never looks in z at all.
    ///
    /// z is degenerate (makeSlab pins every level's slab to z-index 0), so the
    /// refinement ratio in z MUST be 1, not p->ref_vec[2].
    void interp_from_coarse(amrex::iMultiFab& dst, int lev) const
    {
        const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));
        const amrex::IntVect ng = dst.nGrowVect();
        const amrex::BoxArray& dba = dst.boxArray();

        // Coarse patch: one coarsened box per destination box, on the destination's
        // own DistributionMapping, so the MFIter below indexes both with one mfi.
        // Built from a BoxList rather than BoxArray::coarsen because the slab boxes
        // need not be coarsenable — this mirrors how InterpFromCoarseLevel assembles
        // its crse patch from per-box CoarseBox() calls.
        amrex::BoxList bl;
        for (int i = 0; i < dba.size(); ++i)
        {
            bl.push_back(amrex::coarsen(amrex::grow(dba[i], ng), ratio2d));
        }

        amrex::iMultiFab crse_patch(amrex::BoxArray(std::move(bl)), dst.DistributionMap(), 1, 0);
        crse_patch.setVal(0);
        // Cells outside the coarse domain stay 0; InterpFromCoarseLevel leaves them
        // untouched too, since slice_amrex hands it PhysBCFunctNoOp.
        crse_patch.ParallelCopy(m_unique[lev-1], 0, 0, 1, amrex::IntVect(0), amrex::IntVect(0),
                                p->amrex_geometry[lev-1].periodicity());

        for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi)
        {
            auto d = dst.array(mfi);
            auto c = crse_patch.const_array(mfi);
            amrex::ParallelFor(amrex::grow(mfi.validbox(), ng),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    d(i,j,k) = c(amrex::coarsen(amrex::IntVect(AMREX_D_DECL(i,j,k)), ratio2d));
                });
        }
    }

    void checkInvariants(int lev)
    {
        amrex::iMultiFab count(m_unique[lev].boxArray(), m_unique[lev].DistributionMap(), 1, 0);  count.setVal(0);
        count.ParallelAdd(m_owner[lev]);

        // Reduce the max owner-count on the device, then report on the host:
        // amrex::Print()/Abort() cannot be called from inside a GPU kernel.
        amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
        amrex::ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(count); mfi.isValid(); ++mfi)
        {
            auto num  = count.const_array(mfi);
            reduce_op.eval(mfi.validbox(), reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    return { num(i,j,k) };
                });
        }

        int max_count = amrex::get<0>(reduce_data.value());
        amrex::ParallelDescriptor::ReduceIntMax(max_count);

        if(max_count > 1)
        {
            amrex::Print() << "Sliceint invariants violated at lev " << lev << " with max owner count = " << max_count << std::endl;
            amrex::Abort("Sliceint invariants violated.");
        }
    }

    /// Const so the const and non-const accessors share one cache; the cache
    /// state is mutable because filling it is logically const.
    AMREX_FORCE_INLINE void refresh_cache_if_needed() const
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
        {
            // array() on a const FabArray only yields Array4<const int>, so the
            // iMultiFab is un-consted here. Writes through the cache are only
            // reachable via the non-const operator, which requires a non-const
            // sliceint_amrex.
            m_cached_arr4    = const_cast<amrex::iMultiFab&>(m_view[cur_lev])
                                   .atLocalIdx(p->amr_local_fab_idx).array();
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_mfi_idx = cur_idx;
            m_cached_level   = cur_lev;
            m_cached_til_idx = cur_tile_index;
        }
        else if(cur_tile_index != m_cached_til_idx)
        {
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_til_idx = cur_tile_index;
        }
    }

    mutable amrex::Array4<int> m_cached_arr4 = {};
    mutable int m_cached_ox      = 0;
    mutable int m_cached_oy      = 0;
    mutable int m_cached_mfi_idx = -1;
    mutable int m_cached_level   = -1;
    mutable int m_cached_til_idx = -1;

    amrex::Vector<amrex::iMultiFab> m_view;
    amrex::Vector<amrex::iMultiFab> m_unique;
    amrex::Vector<amrex::iMultiFab> m_ksurf;
    amrex::Vector<amrex::iMultiFab> m_owner;
    amrex::Vector<amrex::BoxArray>* m_ba;
    int nlevs;
    lexer *p;
    amrex::IntVect m_ghost = amrex::IntVect(AMREX_D_DECL(0, 0, 0));
};

#endif
#endif
