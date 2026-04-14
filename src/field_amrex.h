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
#ifndef FIELD_AMREX_H_
#define FIELD_AMREX_H_

#include "field.h"
#include "lexer.h"
#include "amrex_bc_func.h"
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_BCUtil.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <initializer_list>
#include <utility>
#include <vector>

#if USE_AMREX
// Create and zero-initialise the combined MultiFab for all fdm fields.
static amrex::Vector<amrex::MultiFab> make_mf(lexer* p, int ncomp)
{
    amrex::Vector<amrex::MultiFab> result(p->nlevs);
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        result[lev].define(p->amrex_box_array[lev],
                           p->amrex_distribution_mapping[lev],
                           ncomp, p->margin);
        result[lev].setVal(0, 0, ncomp, p->margin);
    }
    return result;
}
#endif

class field_amrex : public field
{
public:
    virtual ~field_amrex() = default;

    /*!
     * @copydoc field_base::operator()(int, int, int)
     */
    inline double& operator()(int ii, int jj, int kk) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(ii + m_cached_ox, jj + m_cached_oy, kk + m_cached_oz, 0);
    }

    /*!
     * @copydoc field_base::operator()(int, int, int) const
     */
    inline const double& operator()(int ii, int jj, int kk) const noexcept override final
    {
        refresh_const_cache_if_needed();
        return m_cached_const_arr4(ii + m_cached_const_ox, jj + m_cached_const_oy, kk + m_cached_const_oz, 0);
    }

    /*!
     * @copydoc field_base::operator()(amrex::IntVect, int)
     */
    inline double& operator()(const amrex::IntVect& iv, int comp = 0) noexcept override final
    {
        refresh_cache_if_needed();
        return m_cached_arr4(iv, comp);
    }

    /*!
     * @copydoc field_base::operator()(amrex::IntVect, int) const
     */
    inline const double& operator()(const amrex::IntVect& iv, int comp = 0) const noexcept override final
    {
        refresh_const_cache_if_needed();
        return m_cached_const_arr4(iv, comp);
    }

    /*!
     * @copydoc field_base::setVal()
     */
    void setVal(double val, bool includeGhost = false) override final;

    /*!
     * @copydoc field_base::FillBoundary()
     */
    void FillBoundary() override;

    void FillDomainBoundaryValue(double value, int dir, bool high) override;

    inline amrex::MultiFab& GetMultiFab() noexcept override final {return get_mf(p->level);};
    inline const amrex::MultiFab& GetMultiFab() const noexcept override final {return get_mf_const(p->level);};
    inline amrex::MultiFab& GetMultiFab(int level) noexcept override final {return get_mf(level);};
    inline const amrex::MultiFab& GetMultiFab(int level) const noexcept override final {return get_mf_const(level);};

    /// Returns the shared MultiFab vector pointer (non-null in view mode only).
    amrex::Vector<amrex::MultiFab>* get_shared_mf_vec() noexcept { return m_shared_mf; }


    /// Returns the stagger type of this field.
    amrex_bc_func::DataLocation dataLocation() const
    { return const_params.data_location; }

    // -----------------------------------------------------------------
    // Per-field BCRec update (host side only — no FillPatch).
    // Must be overridden by field1/2/3/4 to use their specific BCDecision.
    // -----------------------------------------------------------------
    virtual void UpdateBCRecs(int gcv) = 0;

    // -----------------------------------------------------------------
    // Batch operations — fill ghost cells for multiple fields that share
    // the same underlying MultiFab in a single AMReX call.
    // -----------------------------------------------------------------

    /// Fill ghost cells (MPI exchange) for a component range of @p shared_mf
    /// in one FillBoundary call instead of one per field.
    static void FillBoundaryBatch(lexer* p,
                                  amrex::Vector<amrex::MultiFab>& shared_mf,
                                  int scomp, int ncomp);

    /// Fill domain boundary conditions for multiple fields that share
    /// @p shared_mf in a single FillPatchSingleLevel / FillPatchTwoLevels call.
    /// Each field must be in view mode (constructed with the view constructor)
    /// and must point to @p shared_mf.
    /// @p scomp is the first component in @p shared_mf belonging to this batch.
    /// Each element of @p fields_and_gcvs is a {field_ptr, gcv} pair so that
    /// each field can be filled with its own ghost-cell variable index.
    static void FillDomainBoundaryBatch(lexer* p,
                                        amrex::Vector<amrex::MultiFab>& shared_mf,
                                        int scomp,
                                        std::initializer_list<std::pair<field_amrex*, int>> fields_and_gcvs,
                                        amrex::Gpu::DeviceVector<amrex::BCRec>& d_bcrec_cache);

protected:
    /// Owning constructor: the field allocates and owns its own MultiFab storage.
    field_amrex(lexer* p, amrex_bc_func::DataLocation data_location);

    /// View constructor: the field is a non-owning view into @p shared_mf at
    /// component @p comp.  The caller must ensure @p shared_mf outlives this object.
    field_amrex(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp,
                amrex_bc_func::DataLocation data_location);

    lexer *p = nullptr;
    amrex::Vector<amrex::MultiFab> mf = {};          ///< owned storage (empty in view mode)
    amrex::Vector<amrex::Vector<amrex::BCRec>> BCRecs = {};

    amrex::Vector<amrex::MultiFab>* m_shared_mf = nullptr; ///< non-owning ptr (view mode only)
    amrex::Vector<amrex::MultiFab> m_alias = {}; ///< 1-component aliases for GetMultiFab() and get_alias() in view mode; empty in owning mode

    // Array4 cache — refreshed once per tile; amortises Array4 construction and
    // tilebox lbound lookup across all cell accesses within the same tile.
    amrex::Array4<amrex::Real> m_cached_arr4 = {};
    int m_cached_ox      = 0;
    int m_cached_oy      = 0;
    int m_cached_oz      = 0;
    int m_cached_mfi_idx = -1;
    int m_cached_level   = -1;
    int m_cached_til_idx = -1;

    mutable amrex::Array4<const amrex::Real> m_cached_const_arr4 = {};
    mutable int m_cached_const_ox      = 0;
    mutable int m_cached_const_oy      = 0;
    mutable int m_cached_const_oz      = 0;
    mutable int m_cached_const_mfi_idx = -1;
    mutable int m_cached_const_level   = -1;
    mutable int m_cached_const_til_idx = -1;

    int m_last_gcv = -1;  ///< gcv from last UpdateBCRecs call; -1 means BCRecs are stale

    /// Cached face BC labels — derived from BCRecs[0][0] and updated only when gcv changes.
    /// Eliminates 6 BCRec reads per FillDomainBoundaryImpl call on the hot path.
    amrex::Array<int, 6> m_cached_face_labels = {};

    /// Device-accessible copy of BCRecs[0], refreshed whenever gcv changes.
    /// Used by the direct-slab ParallelFor path at level 0 to avoid the
    /// PhysBCFunct/GpuBndryFuncFab dispatch layer.
    amrex::Gpu::DeviceVector<amrex::BCRec> m_d_bcrec_lev0 = {};

    template <typename BCDecision>
    void FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision);

    /// Updates BCRecs for all levels using @p bc_decision (host side only).
    /// Called from FillDomainBoundaryImpl and from UpdateBCRecs overrides.
    template <typename BCDecision>
    void UpdateBCRecsImpl(int gcv, const BCDecision& bc_decision);

private:

    /// Refreshes m_cached_arr4 and the tile lbound offset when the current tile or
    /// level has changed.  Called once per operator() invocation; the branch is
    /// almost-always-not-taken (one FAB change per tile loop iteration).
    /// Tile bounds are pre-computed by TILE_LOOP via set_tile_mfi — no
    /// tilebox() or lbound() calls needed here.
    AMREX_FORCE_INLINE void refresh_cache_if_needed()
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_level || cur_idx != m_cached_mfi_idx)
        {
            m_cached_arr4    = get_array(cur_lev,*(p->amr_cell_mfi));
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_oz      = p->amr_tile_lo.z;
            m_cached_mfi_idx = cur_idx;
            m_cached_level   = cur_lev;
            m_cached_til_idx = cur_tile_index;
        }

        if(cur_tile_index != m_cached_til_idx)
        {
            m_cached_ox      = p->amr_tile_lo.x;
            m_cached_oy      = p->amr_tile_lo.y;
            m_cached_oz      = p->amr_tile_lo.z;
            m_cached_til_idx = cur_tile_index;
        }
    }

    /// Refreshes m_cached_const_arr4 when the current tile or level changes.
    AMREX_FORCE_INLINE void refresh_const_cache_if_needed() const
    {
        const int cur_lev = p->level;
        const int cur_idx = p->amr_fab_mfi_idx;
        const int cur_tile_index = p->amr_local_tile_idx;
        if (cur_lev != m_cached_const_level || cur_idx != m_cached_const_mfi_idx)
        {
            m_cached_const_arr4    = get_array_const(cur_lev,*(p->amr_cell_mfi));
            m_cached_const_ox      = p->amr_tile_lo.x;
            m_cached_const_oy      = p->amr_tile_lo.y;
            m_cached_const_oz      = p->amr_tile_lo.z;
            m_cached_const_mfi_idx = cur_idx;
            m_cached_const_level   = cur_lev;
            m_cached_const_til_idx = cur_tile_index;
        }

        if(cur_tile_index != m_cached_const_til_idx)
        {
            m_cached_const_ox      = p->amr_tile_lo.x;
            m_cached_const_oy      = p->amr_tile_lo.y;
            m_cached_const_oz      = p->amr_tile_lo.z;
            m_cached_const_til_idx = cur_tile_index;
        }
    }

    /// Returns the 1-component alias for @p level (alias in view mode, owned otherwise).
    AMREX_FORCE_INLINE amrex::MultiFab& get_mf(int level) noexcept
    { return m_shared_mf ? m_alias[level] : mf[level]; }

    /// Const overload — Returns the 1-component alias for @p level (alias in view mode, owned otherwise).
    AMREX_FORCE_INLINE const amrex::MultiFab& get_mf_const(int level) const noexcept
    { return m_shared_mf ? m_alias[level] : mf[level]; }

    /// Returns the storage MultiFab for @p level (shared in view mode, owned otherwise).
    AMREX_FORCE_INLINE amrex::MultiFab& get_mf_shared(int level) noexcept
    { return m_shared_mf ? (*m_shared_mf)[level] : mf[level]; }

    /// Const overload — used by the const cache refresh to access data without mutation.
    AMREX_FORCE_INLINE const amrex::MultiFab& get_mf_const_shared(int level) const noexcept
    { return m_shared_mf ? (*m_shared_mf)[level] : mf[level]; }

    /// Returns the Array4 for the current tile of @p level.
    AMREX_FORCE_INLINE amrex::Array4<amrex::Real> get_array(int level, amrex::MFIter& mfi) noexcept
    { return get_mf(level).array(mfi); }

    /// Const overload — used by the const cache refresh to access data without mutation.
    AMREX_FORCE_INLINE const amrex::Array4<const amrex::Real> get_array_const(int level, amrex::MFIter& mfi) const noexcept
    { return get_mf_const(level).const_array(mfi); }

    /// Shifts face data inward at the high-end boundary for face-staggered fields.
    static void ShiftBigBoundaryFaceInward(amrex::MultiFab& mf_in,
                                           amrex_bc_func::DataLocation data_location,
                                           const amrex::Geometry& geom);

    /// Fills ghost-cell slabs on all 6 domain faces via direct ParallelFor calls.
    /// Shared by FillDomainBoundaryImpl (single-component, scomp=0) and
    /// FillDomainBoundaryBatch (multi-component, scomp=batch offset).
    template<typename BCDecision>
    static void fill_boundary_slabs(
        amrex::MultiFab& mf_lev,
        int scomp,
        int ncomp,
        const amrex::BCRec* bcrec_d,
        const amrex_bc_func::MyExtBCFillField<BCDecision>& fill,
        const amrex::GeometryData& geom_data,
        const amrex::Box& dom);

    /// Build a flat BCRec vector for @p lev by taking the first BCRec from
    /// each field in @p fields, in order.
    static amrex::Vector<amrex::BCRec>
    make_combined_bcrecs(int lev,
                         std::initializer_list<std::pair<field_amrex*, int>> fields_and_gcvs)
    {
        amrex::Vector<amrex::BCRec> result;
        result.reserve(fields_and_gcvs.size());
        for (const auto& [f, gcv] : fields_and_gcvs)
            if (!f->BCRecs[lev].empty())
                result.push_back(f->BCRecs[lev][0]);
        return result;
    }

    const amrex_bc_func::ConstMyExtBCFillFieldParams const_params = {};
};

// =========================================================================
// Namespace with BC decision helper utilities (used by field1-4 constructors)
// =========================================================================
namespace field_amrex_detail
{
    using Label = amrex_bc_func::BoundaryConditionTypeLabel;

    enum class OutflowAxis { U, V, W };

    inline Label base_parallel_label(int B20)
    {
        switch (B20)
        {
            case 1: return Label::NEUMANN;
            case 2: case 4: return Label::NOSLIP;
            case 3: return Label::DIRICHLET_ORTH;
            default: return Label::NEUMANN;
        }
    }

    inline std::pair<Label, Label> compute_parallel_u_v(lexer* p)
    {
        Label u = base_parallel_label(p->B20);
        Label v = base_parallel_label(p->B20);
        if (p->B23 == 2)
        {
            u = Label::DIRICHLET_PARA_REFLECT;
            v = Label::DIRICHLET_PARA_REFLECT;
        }

        if (p->A217 == 1 && p->A10 == 2)
        {
            u = Label::NEUMANN;
            v = Label::NEUMANN;
        }
        else if (p->A217 == 2 && p->A10 == 2)
        {
            u = Label::NOSLIP;
            v = Label::NOSLIP;
        }

        return {u, v};
    }

    inline Label compute_parallel_w(lexer* p)
    {
        if (p->B23 == 2)
            return Label::DIRICHLET_PARA_REFLECT;
        return base_parallel_label(p->B20);
    }

    inline Label compute_orth_label(lexer* p)
    {
        return (p->B23 == 2) ? Label::DIRICHLET_ORTH_REFLECT : Label::DIRICHLET_ORTH;
    }

    inline Label compute_inflow_label(lexer* p)
    {
        if (p->I230 > 0 || p->B98 >= 3 || p->B60 > 0)
            return Label::NONE;
        return Label::DIRICHLET_ORTH;
    }

    inline Label compute_outflow_label(lexer* p, OutflowAxis axis)
    {
        switch (p->B75)
        {
            case 1:
                return Label::NEUMANN;
            case 2:
                return Label::OUTFLOWBC;
            case 3:
                return (axis == OutflowAxis::U) ? Label::NONE : Label::OUTFLOWBC;
            default:
                return Label::NEUMANN;
        }
    }

    inline bool compute_awa_label(lexer* p)
    {
        return p->B99 >= 3;
    }

    inline bool compute_gclabel_outflow(lexer* p)
    {
        return !(p->B60 == 3 || p->B60 == 4);
    }

    inline bool compute_i10_enabled(lexer* p)
    {
        return p->I10 == 1;
    }
}

// =========================================================================
// UpdateBCRecsImpl — host-side BCRec population (defined inline in header
// because it is templated on BCDecision)
// =========================================================================
template <typename BCDecision>
void field_amrex::UpdateBCRecsImpl(int gcv, const BCDecision& bc_decision)
{
    // Direction integer codes used by the BC decision infrastructure
    const int x_neg = 1;  // Dir::X_NEG
    const int x_pos = 4;  // Dir::X_POS
    const int y_neg = 3;  // Dir::Y_NEG
    const int y_pos = 2;  // Dir::Y_POS
    const int z_neg = 5;  // Dir::Z_NEG
    const int z_pos = 6;  // Dir::Z_POS

    LEVEL_LOOP
    {
        if (BCRecs[p->level].empty()) continue;

        auto& bc = BCRecs[p->level][0];

        bc.setLo(0, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[0], x_neg)));
        bc.setHi(0, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[1], x_pos)));
        bc.setLo(1, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[2], y_neg)));
        bc.setHi(1, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[3], y_pos)));
        bc.setLo(2, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[4], z_neg)));
        bc.setHi(2, static_cast<int>(bc_decision.evaluate(gcv, const_params.bc_values[5], z_pos)));
    }
    m_last_gcv = gcv;
}

// =========================================================================
// fill_boundary_slabs — fills the 6 ghost-cell face slabs via ParallelFor.
// Extracted from FillDomainBoundaryImpl and FillDomainBoundaryBatch to
// eliminate the duplicated per-face slab logic.
// =========================================================================
template<typename BCDecision>
void field_amrex::fill_boundary_slabs(
    amrex::MultiFab& mf_lev,
    int scomp,
    int ncomp,
    const amrex::BCRec* bcrec_d,
    const amrex_bc_func::MyExtBCFillField<BCDecision>& fill,
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
            const amrex::Box slab(
                amrex::IntVect(fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)),
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
                amrex::IntVect(fabbox.bigEnd(0), fabbox.bigEnd(1), fabbox.bigEnd(2)));
            amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                     amrex::Real(0), bcrec_d, 0, 0);
            });
        }
        // Z-lo slab
        if (fabbox.smallEnd(2) < dom.smallEnd(2))
        {
            const amrex::Box slab(
                amrex::IntVect(fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)),
                amrex::IntVect(fabbox.bigEnd(0), fabbox.bigEnd(1), dom.smallEnd(2)-1));
            amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                     amrex::Real(0), bcrec_d, 0, 0);
            });
        }
        // Z-hi slab
        if (fabbox.bigEnd(2) > dom.bigEnd(2))
        {
            const amrex::Box slab(
                amrex::IntVect(fabbox.smallEnd(0), fabbox.smallEnd(1), dom.bigEnd(2)+1),
                fabbox.bigEnd());
            amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                fill(amrex::IntVect(i,j,k), arr, scomp, ncomp, geom_data,
                     amrex::Real(0), bcrec_d, 0, 0);
            });
        }
    }
}

// =========================================================================
// FillDomainBoundaryImpl — direct domain BC fill (no MPI exchange here)
//
// Level 0: iterate only the 6 ghost-cell slabs (X+/-, Y+/-, Z+/-) and call
//   MyExtBCFillField directly — no PhysBCFunct/GpuBndryFuncFab overhead,
//   no double-fill, no interior-cell processing.  MPI exchange is done
//   separately by FillBoundary() / start*MPI calls.
//
// Level > 0: unchanged FillPatchTwoLevels path (handles MPI + interpolation).
// =========================================================================
template <typename BCDecision>
void field_amrex::FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision)
{
    if (gcv != m_last_gcv)
    {
        UpdateBCRecsImpl(gcv, bc_decision);
        // Cache face labels — BCRecs[0][0] is the same for all levels; only
        // changes when gcv changes, so avoid re-reading on every call.
        if (!BCRecs.empty() && !BCRecs[0].empty())
        {
            const auto& bc = BCRecs[0][0];
            m_cached_face_labels[0] = bc.lo(2);  // Z_NEG
            m_cached_face_labels[1] = bc.hi(2);  // Z_POS
            m_cached_face_labels[2] = bc.lo(0);  // X_NEG
            m_cached_face_labels[3] = bc.hi(0);  // X_POS
            m_cached_face_labels[4] = bc.lo(1);  // Y_NEG
            m_cached_face_labels[5] = bc.hi(1);  // Y_POS

            // Refresh device BCRec copy for the direct-slab ParallelFor path.
            m_d_bcrec_lev0.resize(BCRecs[0].size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             BCRecs[0].begin(), BCRecs[0].end(),
                             m_d_bcrec_lev0.begin());
        }
    }

    // Build functor on the stack — cheap (~100 B copy).
    // Ui/Uo/dt are time-varying and always read fresh from lexer.
    amrex_bc_func::MyExtBCFillFieldParams params(p->Ui, p->Uo, p->dt, m_cached_face_labels);

    const auto fill = amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params};
    // bf (GpuBndryFuncFab) is only needed for level > 0 (FillPatchTwoLevels path).
    // For single-level runs it is never used, so construction is deferred into the
    // else branch below to avoid the overhead on every start1/2/3/4 call.

    LEVEL_LOOP
    {
        auto& mf_lev = get_mf(p->level);
        if(p->level==0)
        {
            mf_lev.FillBoundary(0, 1, p->amrex_geometry[p->level].periodicity());

            // Direct slab ParallelFor — bypasses PhysBCFunct/GpuBndryFuncFab dispatch.
            // For each FAB that touches a domain boundary, launch one ParallelFor per
            // face direction covering only that face's ghost-cell slab.  Edge/corner
            // cells lie in multiple slabs and are written more than once; the functor is
            // deterministic so the repeated writes are idempotent.
            const amrex::GeometryData geom_data = p->amrex_geometry[p->level].data();
            const amrex::Box          dom       = p->amrex_geometry[p->level].Domain();

            fill_boundary_slabs(mf_lev, 0, mf_lev.nComp(), m_d_bcrec_lev0.data(),
                                fill, geom_data, dom);
        }
        else
        {
            auto& mf_coarse = get_mf(p->level-1);
            // Multi-level: FillPatchTwoLevels handles MPI exchange and
            // coarse-to-fine interpolation; keep the existing PhysBCFunct path.
            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> bf(fill);

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> cphysbcf(
                p->amrex_geometry[p->level-1], BCRecs[p->level-1], bf);

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> fphysbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], bf);

            amrex::Interpolater* mapper = &amrex::cell_cons_interp;
            const amrex::IntVect ref_vec = p->ref_vec;

            amrex::FillPatchTwoLevels(mf_lev, amrex::Real(p->simtime),
                                        {&(mf_coarse)}, {amrex::Real(p->simtime)},
                                        {&(mf_lev)}, {amrex::Real(p->simtime)},
                                        0, 0, mf_lev.nComp(), p->amrex_geometry[p->level-1], p->amrex_geometry[p->level],
                                        cphysbcf, 0,
                                        fphysbcf, 0,
                                        ref_vec,
                                        mapper,
                                        BCRecs[p->level], 0);
        }

        ShiftBigBoundaryFaceInward(mf_lev, const_params.data_location, p->amrex_geometry[p->level]);
    }
}

// =========================================================================
// FillDomainBoundaryBatch — one FillPatch for all fields in the batch.
// The BCDecision type is used only to satisfy the MyExtBCFillField template;
// actual per-component BC behaviour is encoded in the BCRecs populated by
// each field's UpdateBCRecs() call, so any BCDecision type is equivalent.
// Use Field4BcDecision as a concrete, always-available phantom type.
// =========================================================================
inline void field_amrex::FillDomainBoundaryBatch(
    lexer* p,
    amrex::Vector<amrex::MultiFab>& shared_mf,
    int scomp,
    std::initializer_list<std::pair<field_amrex*, int>> fields_and_gcvs,
    amrex::Gpu::DeviceVector<amrex::BCRec>& d_bcrec_cache)
{
    // Step 1 — populate each field's BCRecs on the host using its own gcv
    for (auto& [f, gcv] : fields_and_gcvs)
        if (f->m_last_gcv != gcv)
            f->UpdateBCRecs(gcv);

    const int ncomp = static_cast<int>(fields_and_gcvs.size());

    // Step 2 — build a const_params representative for this group
    // (bc_values, heat_values, y_dimension, data_location are the same for every
    //  field in a stagger group since they all come from the same lexer p)
    amrex_bc_func::ConstMyExtBCFillFieldParams const_params{
        {p->bcside1, p->bcside4, p->bcside3, p->bcside2, p->bcside5, p->bcside6},
        {p->H61_T, p->H64_T, p->H63_T, p->H62_T, p->H65_T, p->H66_T},
        bool(p->j_dir), amrex_bc_func::DataLocation::CELL_CENTERED};

    amrex_bc_func::MyExtBCFillFieldParams params;
    params.Ui = p->Ui;
    params.Uo = p->Uo;
    params.dt = p->dt;

    // Step 3 — one FillPatch call per level covering [scomp, scomp+ncomp)
    using PhantomDecision = amrex_bc_func::Field4BcDecision;

    LEVEL_LOOP
    {
        // Build combined BCRecs (ncomp entries) from each field's pre-populated BCRecs
        auto combined      = make_combined_bcrecs(p->level, fields_and_gcvs);

        if (p->level == 0)
        {
            // face_labels is left default: ncomp > 1 falls through to BCRec pointer reads.
            const auto fill = amrex_bc_func::MyExtBCFillField<PhantomDecision>{const_params, params};

            const amrex::GeometryData geom_data = p->amrex_geometry[0].data();
            const amrex::Box          dom       = p->amrex_geometry[0].Domain();

            // Reuse the caller-owned device buffer — resize only when the batch
            // size changes (no-op for the common steady-state case).
            d_bcrec_cache.resize(ncomp);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             combined.begin(), combined.end(),
                             d_bcrec_cache.begin());
            const amrex::BCRec* bcrec_d = d_bcrec_cache.data();

            auto& mf_lev = shared_mf[p->level];
            fill_boundary_slabs(mf_lev, scomp, ncomp, bcrec_d, fill, geom_data, dom);
            amrex::Gpu::streamSynchronize();
        }
        else
        {
            auto combined_prev = make_combined_bcrecs(p->level - 1, fields_and_gcvs);

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<PhantomDecision>> cbf(
                amrex_bc_func::MyExtBCFillField<PhantomDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<PhantomDecision>>> cphysbcf(
                p->amrex_geometry[p->level-1], combined_prev, cbf);

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<PhantomDecision>> fbf(
                amrex_bc_func::MyExtBCFillField<PhantomDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<PhantomDecision>>> fphysbcf(
                p->amrex_geometry[p->level], combined, fbf);

            amrex::Interpolater* mapper = &amrex::cell_cons_interp;
            const amrex::IntVect ref_vec = p->ref_vec;

            amrex::FillPatchTwoLevels(shared_mf[p->level], amrex::Real(p->simtime),
                                      {&shared_mf[p->level-1]}, {amrex::Real(p->simtime)},
                                      {&shared_mf[p->level]},   {amrex::Real(p->simtime)},
                                      scomp, scomp, ncomp,
                                      p->amrex_geometry[p->level-1], p->amrex_geometry[p->level],
                                      cphysbcf, 0,
                                      fphysbcf, 0,
                                      ref_vec, mapper,
                                      combined, 0);
        }

        // Shift face data at the high-end boundary per-field using each
        // field's own stagger type and its 1-component alias.
        for (auto& [f, gcv] : fields_and_gcvs)
            ShiftBigBoundaryFaceInward(f->GetMultiFab(),
                                       f->dataLocation(),
                                       p->amrex_geometry[p->level]);
    }
}

#endif
#endif
