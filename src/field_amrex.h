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
#include <AMReX_Vector.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <initializer_list>
#include <utility>
#include <vector>

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
        else if(cur_tile_index != m_cached_til_idx)
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
        else if(cur_tile_index != m_cached_const_til_idx)
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

    void ShiftBigBoundaryFaceInward(amrex::MultiFab& mf_in, int p_level)
    {
        int dir = -1;
        if (const_params.data_location == amrex_bc_func::DataLocation::FACE_X) dir = 0;
        else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Y) dir = 1;
        else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Z) dir = 2;

        if (dir == -1) return;

        const auto& geom = p->amrex_geometry[p_level];
        const int domain_hi = geom.Domain().bigEnd(dir);

        for (amrex::MFIter mfi(mf_in); mfi.isValid(); ++mfi)
        {
            const amrex::Box& valid_box = mfi.validbox();

            // Check if this box touches the high boundary in the specific direction
            if (valid_box.bigEnd(dir) == domain_hi)
            {
                const amrex::Box& box = mfi.fabbox();
                amrex::Array4<amrex::Real> const& arr = mf_in.array(mfi);

                int start = domain_hi;
                int end = box.bigEnd(dir) - 1;

                // Define a box collapsed to the start plane for iteration
                amrex::Box para_box = box;
                para_box.setSmall(dir, start);
                para_box.setBig(dir, start);

                if (dir == 0)
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int /*i dummy*/, int j, int k, int n)
                    {
                        for (int i = start; i <= end; ++i)
                        {
                            arr(i, j, k, n) = arr(i + 1, j, k, n);
                        }
                    });
                }
                else if (dir == 1)
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int i, int /*j dummy*/, int k, int n)
                    {
                        for (int j = start; j <= end; ++j)
                        {
                            arr(i, j, k, n) = arr(i, j + 1, k, n);
                        }
                    });
                }
                else // dir == 2
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k dummy*/, int n)
                    {
                        for (int k = start; k <= end; ++k)
                        {
                            arr(i, j, k, n) = arr(i, j, k + 1, n);
                        }
                    });
                }
            }
        }
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
        }
    }

    // Build functor on the stack — cheap (~100 B copy).
    // Ui/Uo/dt are time-varying and always read fresh from lexer.
    amrex_bc_func::MyExtBCFillFieldParams params(p->Ui, p->Uo, p->dt, m_cached_face_labels);

    LEVEL_LOOP
    {
        auto& mf_lev = get_mf(p->level);
        if(p->level==0)
        {
            mf_lev.FillBoundary(0, 1, p->amrex_geometry[p->level].periodicity());

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> bf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> physbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], bf);

            physbcf(mf_lev, 0, mf_lev.nComp(), mf_lev.nGrowVect(), amrex::Real(p->simtime), 0);
        }
        else
        {
            auto& mf_coarse = get_mf(p->level-1);
            // Multi-level: FillPatchTwoLevels handles MPI exchange and
            // coarse-to-fine interpolation; keep the existing PhysBCFunct path.
            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> cbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> cphysbcf(
                p->amrex_geometry[p->level-1], BCRecs[p->level-1], cbf);

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> fbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> fphysbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], fbf);

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

        ShiftBigBoundaryFaceInward(mf_lev, p->level);
    }
}

#endif
#endif
