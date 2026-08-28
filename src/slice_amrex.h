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

#if USE_AMREX
#ifndef SLICE_AMREX_H_
#define SLICE_AMREX_H_

#include "slice_amrex_prim.h"
#include "lexer.h"
#include "amrex_bc_func2D.h"
#include <AMReX_Reduce.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>

class slice_amrex : public slice_amrex_prim
{
public:
    slice_amrex(lexer *pp, DataLocation data_location)
        : const_params({pp->bcside1, pp->bcside4, pp->bcside3, pp->bcside2},
                       pp->j_dir, data_location, sqrt(pp->W20*pp->W20+pp->W21*pp->W21+pp->W22*pp->W22))
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

        // One BCRec per level; the slice carries a single component.
        BCRecs.resize(nlevs);
        for (auto& bc_rec : BCRecs)
            bc_rec.resize(1);

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
    virtual ~slice_amrex()
    {
        if(p)
        {
            p->deregister_slice(this);
        }
    };

    slice_amrex(const slice_amrex&) = delete;
    slice_amrex& operator=(const slice_amrex&) = delete;
    slice_amrex(slice_amrex&&) = delete;
    slice_amrex& operator=(slice_amrex&&) = delete;

    inline double& operator()(int ii, int jj) noexcept
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, 0, 0);
    };

    inline const double& operator()(int ii, int jj) const noexcept
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, 0, 0);
    };

    void regrid() override final
    {
        const int old_nlevs = nlevs;
        nlevs = p->nlevs;
        m_view.resize(nlevs);
        m_unique.resize(nlevs);
        m_ksurf.resize(nlevs);
        m_owner.resize(nlevs);

        // Keep BCRecs level-aligned; new levels start with one (stale) BCRec, and
        // m_last_gcv is invalidated below so the next fill repopulates them.
        BCRecs.resize(nlevs);
        for (int lev = old_nlevs; lev < nlevs; ++lev)
            BCRecs[lev].resize(1);
        if (nlevs != old_nlevs) m_last_gcv = -1;

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
                amrex::MultiFab new_unique(ba2d_unique, dm2d, 1, m_ghost);
                new_unique.setVal(0, 0, 1, m_ghost);

                if (lev > 0)
                {
                    // z degenerate -> ratio 1 in z; pc_interp for the same reason as
                    // everywhere else here (single-cell stencil, never looks in z).
                    const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));
                    amrex::PhysBCFunctNoOp nobc;
                    amrex::Vector<amrex::BCRec> bcr(1);
                    amrex::InterpFromCoarseLevel(new_unique, amrex::Real(p->simtime),
                        m_unique[lev-1], 0, 0, 1,
                        p->amrex_geometry[lev-1], p->amrex_geometry[lev],
                        nobc, 0, nobc, 0,
                        ratio2d, &amrex::pc_interp, bcr, 0);
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

    // -----------------------------------------------------------------
    // Per-slice BCRec update (host side only — no FillPatch).
    // Must be overridden by the concrete slice types to supply their BCDecision.
    // -----------------------------------------------------------------
    virtual void UpdateBCRecs(int gcv) = 0;
    virtual void FillDomainBoundary(int gcv) = 0;

protected:
    const amrex_bc_func2D::ConstMyExtBCFillSliceParams const_params = {};

    // =====================================================================
    // UpdateBCRecsImpl — 2D analogue of field_amrex::UpdateBCRecsImpl.
    // Only x/y have physical boundaries; the slab is degenerate in z.
    // =====================================================================
    template <typename BCDecision>
    void UpdateBCRecsImpl(int gcv, const BCDecision& bc_decision)
    {
        // Direction integer codes used by the BC decision infrastructure
        // (amrex_bc_func2D::Dir is private, so the raw codes are repeated here
        // exactly as field_amrex::UpdateBCRecsImpl does).
        const int x_neg = 1;  // Dir::X_NEG
        const int x_pos = 4;  // Dir::X_POS
        const int y_neg = 3;  // Dir::Y_NEG
        const int y_pos = 2;  // Dir::Y_POS

        for (int lev = 0; lev < nlevs; ++lev)
        {
            if (BCRecs[lev].empty()) continue;

            auto& bc = BCRecs[lev][0];

            bc.setLo(0, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[0], x_neg)));
            bc.setHi(0, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[1], x_pos)));
            bc.setLo(1, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[2], y_neg)));
            bc.setHi(1, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[3], y_pos)));
            // z has no physical boundary on a slab, but BCRec is AMREX_SPACEDIM-sized:
            // set it explicitly so the entry is never read uninitialised.
            bc.setLo(2, 0);
            bc.setHi(2, 0);
        }
        m_last_gcv = gcv;
    }

    // =====================================================================
    // fill_boundary_slabs — fills the 4 ghost-cell face slabs via ParallelFor.
    // 2D analogue of field_amrex::fill_boundary_slabs: the Z-lo/Z-hi slabs are
    // gone because the slice carries no z ghost cells (m_ghost[2] == 0), so
    // every fabbox spans exactly z = [0,0].
    // =====================================================================
    template <typename BCDecision>
    static void fill_boundary_slabs(
        amrex::MultiFab& mf_lev,
        int scomp,
        int ncomp,
        const amrex::BCRec* bcrec_d,
        const amrex_bc_func2D::MyExtBCFillSlice<BCDecision>& fill,
        const amrex::GeometryData& geom_data,
        const amrex::Box& dom)
    {
        for (amrex::MFIter mfi(mf_lev); mfi.isValid(); ++mfi)
        {
            const amrex::Box& fabbox = mfi.fabbox();
            if (dom.contains(fabbox)) continue;

            auto arr = mf_lev.array(mfi);

            // X-lo slab
            if (fabbox.smallEnd(0) < dom.smallEnd(0))
            {
                const amrex::Box slab(fabbox.smallEnd(),
                    amrex::IntVect(dom.smallEnd(0)-1, fabbox.bigEnd(1), fabbox.bigEnd(2)));
                amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                         amrex::Real(0), bcrec_d, 0, 0);
                });
            }
            // X-hi slab
            if (fabbox.bigEnd(0) > dom.bigEnd(0))
            {
                const amrex::Box slab(
                    amrex::IntVect(dom.bigEnd(0)+1, fabbox.smallEnd(1), fabbox.smallEnd(2)),
                    fabbox.bigEnd());
                amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                         amrex::Real(0), bcrec_d, 0, 0);
                });
            }
            // Y-lo slab
            if (fabbox.smallEnd(1) < dom.smallEnd(1))
            {
                const amrex::Box slab(fabbox.smallEnd(),
                    amrex::IntVect(fabbox.bigEnd(0), dom.smallEnd(1)-1, fabbox.bigEnd(2)));
                amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                         amrex::Real(0), bcrec_d, 0, 0);
                });
            }
            // Y-hi slab
            if (fabbox.bigEnd(1) > dom.bigEnd(1))
            {
                const amrex::Box slab(
                    amrex::IntVect(fabbox.smallEnd(0), dom.bigEnd(1)+1, fabbox.smallEnd(2)),
                    fabbox.bigEnd());
                amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                         amrex::Real(0), bcrec_d, 0, 0);
                });
            }
        }
    }

    // =====================================================================
    // FillDomainBoundaryImpl — 2D analogue of field_amrex::FillDomainBoundaryImpl.
    //
    // Same lev==0 / lev>0 split as the field version, but the view/unique separation
    // changes what each branch is allowed to do.  m_view is deliberately OVERLAPPING
    // (stacked-in-z 3D boxes flatten onto the same column), so it can never be a
    // FillBoundary/FillPatch *source* — but it is fine as a *destination*, since each
    // FAB is then filled independently from m_unique, a removeOverlap partition.
    //
    //   lev == 0 : makeUnique() broadcasts m_unique's valid cells into m_view's valid
    //              cells and interior ghosts (standing in for field_amrex's
    //              FillBoundary, which is illegal here), then fill_boundary_slabs()
    //              writes the domain-boundary ghosts.  That slab fill is purely local
    //              — each ghost reads its own FAB's interior — so overlapping columns
    //              all compute the same value from interiors makeUnique just made
    //              consistent.
    //   lev >  0 : FillPatchTwoLevels fills m_view from the m_unique partitions:
    //              fine data where this level covers, coarse interpolation for the
    //              coarse-fine ghost ring, physical BCs at the domain boundary.
    //
    // Note fillHoles() is still required and is NOT subsumed by FillPatchTwoLevels:
    // FillPatch fills by COVERAGE, and a hole is a column m_unique[lev] does cover —
    // it merely holds 0 there — so FillPatch would happily propagate the zero.
    // =====================================================================
    template <typename BCDecision>
    void FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision)
    {
        if (gcv != m_last_gcv)
        {
            UpdateBCRecsImpl(gcv, bc_decision);
            // Cache face labels — BCRecs[0][0] is the same for all levels; only
            // changes when gcv changes.
            if (!BCRecs.empty() && !BCRecs[0].empty())
            {
                const auto& bc = BCRecs[0][0];
                m_cached_face_labels[0] = bc.lo(0);  // X_NEG
                m_cached_face_labels[1] = bc.hi(0);  // X_POS
                m_cached_face_labels[2] = bc.lo(1);  // Y_NEG
                m_cached_face_labels[3] = bc.hi(1);  // Y_POS

                // Refresh device BCRec copy for the direct-slab ParallelFor path.
                m_d_bcrec_lev0.resize(BCRecs[0].size());
                amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                                 BCRecs[0].begin(), BCRecs[0].end(),
                                 m_d_bcrec_lev0.begin());
            }
        }

        // Build functor on the stack — Ui/Uo/dt/wd are time-varying, always read fresh.
        const amrex_bc_func2D::MyExtBCFillSliceParams params(p->Ui, p->Uo, p->dt, p->wd,
                                                             m_cached_face_labels);
        const auto fill = amrex_bc_func2D::MyExtBCFillSlice<BCDecision>{const_params, params};

        for (int lev = 0; lev < nlevs; ++lev)
        {
            // Reduce to m_unique, fill un-owned columns from lev-1, broadcast back:
            // gives m_view its valid data plus interior (same-level) x/y ghosts.
            makeUnique(lev);

            if (lev == 0)
            {
                // Direct slab ParallelFor — bypasses PhysBCFunct/GpuBndryFuncFab dispatch.
                // (field_amrex FillBoundary()s its own MultiFab at this point; the slice
                // cannot — m_view is overlapping — so makeUnique's broadcast from the
                // m_unique partition plays that role instead.)
                const amrex::GeometryData geom_data = p->amrex_geometry[0].data();
                const amrex::Box          dom       = p->amrex_geometry[0].Domain();

                fill_boundary_slabs(m_view[0], 0, 1, m_d_bcrec_lev0.data(),
                                    fill, geom_data, dom);
            }
            else
            {
                // Fill the view straight from the unique partitions: fine data where this
                // level covers, coarse interpolation for the coarse-fine ghost ring, and
                // physical BCs at the domain boundary — all in one call.
                //
                // The overlapping view is safe as the DESTINATION because every FAB is
                // filled independently from m_unique, which is a proper (removeOverlap)
                // partition; there is no fine-fine exchange among the view's own boxes
                // for the overlap to make ambiguous.  This supersedes the broadcast
                // makeUnique just did for valid cells and interior ghosts.
                amrex::GpuBndryFuncFab<amrex_bc_func2D::MyExtBCFillSlice<BCDecision>> bf(fill);

                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func2D::MyExtBCFillSlice<BCDecision>>>
                    cphysbcf(p->amrex_geometry[lev-1], BCRecs[lev-1], bf);
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func2D::MyExtBCFillSlice<BCDecision>>>
                    fphysbcf(p->amrex_geometry[lev], BCRecs[lev], bf);

                // z is degenerate — every level's slab sits at z-index 0 — so the
                // refinement ratio in z MUST be 1.
                //
                // pc_interp, not cell_cons_interp, and not only for consistency with
                // fillHoles: a slope-based interpolater grows its coarse patch by one
                // cell in every direction, so it would read coarse cells at z = +-1 that
                // the slab does not have (m_ghost[2] == 0).  Those reads are
                // uninitialised, and although a z-ratio of 1 gives the z-slope zero
                // weight, NaN * 0 is still NaN.  pc_interp's stencil is a single coarse
                // cell, so it never looks in z at all.
                const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));
                amrex::Interpolater* mapper = &amrex::pc_interp;

                amrex::FillPatchTwoLevels(m_view[lev], amrex::Real(p->simtime),
                                          {&(m_unique[lev-1])}, {amrex::Real(p->simtime)},
                                          {&(m_unique[lev])},   {amrex::Real(p->simtime)},
                                          0, 0, 1,
                                          p->amrex_geometry[lev-1], p->amrex_geometry[lev],
                                          cphysbcf, 0,
                                          fphysbcf, 0,
                                          ratio2d,
                                          mapper,
                                          BCRecs[lev], 0);
            }
        }
    }

private:
    void buildOwnerMask(int lev)
    {
        const auto& ba3d = (*m_ba)[lev];
        for (amrex::MFIter mfi(m_view[lev]); mfi.isValid(); ++mfi)
        {
            const amrex::Box& b3 = ba3d[mfi.index()];               // the un-flattened box
            auto own  = m_owner[lev].array(mfi);
            auto ksurf = m_ksurf[lev].array(mfi);
            // ksurf(i,j): surface cell index from phi/eta (interp from coarse right after regrid)
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

        // Interpolate the coarser slice onto this level's unique layout. z is
        // degenerate (makeSlab pins every level's slab to z-index 0), so the
        // refinement ratio in z MUST be 1, not _p->ref_vec[2].
        amrex::MultiFab coarse_fill(m_unique[lev].boxArray(),
                                    m_unique[lev].DistributionMap(), 1, 0);
        const amrex::IntVect ratio2d(AMREX_D_DECL(p->ref_vec[0], p->ref_vec[1], 1));

        amrex::PhysBCFunctNoOp nobc;
        amrex::Interpolater*   mapper = &amrex::pc_interp;   // injection: hole column takes coarse column's value
        amrex::Vector<amrex::BCRec> bcr(1);                 // unused with PhysBCFunctNoOp

        amrex::InterpFromCoarseLevel(coarse_fill, amrex::Real(p->simtime),
            m_unique[lev-1], 0, 0, 1,
            p->amrex_geometry[lev-1], p->amrex_geometry[lev],
            nobc, 0, nobc, 0,
            ratio2d, mapper, bcr, 0);

        // Blend: keep reduced (owned) values; take interpolated coarse value in holes.
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
            amrex::Print() << "Slice invariants violated at lev " << lev << " with max owner count = " << max_count << std::endl;
            amrex::Abort("Slice invariants violated.");
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
            // array() on a const FabArray only yields Array4<const Real>, so the
            // MultiFab is un-consted here. Writes through the cache are only
            // reachable via the non-const operator, which requires a non-const
            // slice_amrex.
            m_cached_arr4    = const_cast<amrex::MultiFab&>(m_view[cur_lev])
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

    mutable amrex::Array4<amrex::Real> m_cached_arr4 = {};
    mutable int m_cached_ox      = 0;
    mutable int m_cached_oy      = 0;
    mutable int m_cached_mfi_idx = -1;
    mutable int m_cached_level   = -1;
    mutable int m_cached_til_idx = -1;

    amrex::Vector<amrex::Vector<amrex::BCRec>> BCRecs = {};

    int m_last_gcv = -1;  ///< gcv from last UpdateBCRecs call; -1 means BCRecs are stale

    /// Cached face BC labels — derived from BCRecs[0][0], updated only when gcv changes.
    /// Indexed as: 0=X_NEG  1=X_POS  2=Y_NEG  3=Y_POS
    amrex::Array<int,4> m_cached_face_labels = {};

    /// Device-accessible copy of BCRecs[0], refreshed whenever gcv changes.
    amrex::Gpu::DeviceVector<amrex::BCRec> m_d_bcrec_lev0 = {};

    amrex::Vector<amrex::MultiFab> m_view;
    amrex::Vector<amrex::MultiFab> m_unique;
    amrex::Vector<amrex::MultiFab> m_ksurf;
    amrex::Vector<amrex::iMultiFab> m_owner;
    amrex::Vector<amrex::BoxArray>* m_ba;
    int nlevs;
    lexer *p;
    amrex::IntVect m_ghost = amrex::IntVect(AMREX_D_DECL(0, 0, 0));
};

#endif
#endif
