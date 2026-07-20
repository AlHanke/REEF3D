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
#include "field_amrex.h"
#include "lexer.h"
#include "amrex_bc_func.h"
#include <AMReX_BCUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>

// ---------------------------------------------------------------------------
// Owning constructor
// ---------------------------------------------------------------------------
field_amrex::field_amrex(lexer* p, amrex_bc_func::DataLocation data_location)
    : const_params({p->bcside1, p->bcside4, p->bcside3, p->bcside2, p->bcside5, p->bcside6},
                   {p->H61_T, p->H64_T, p->H63_T, p->H62_T, p->H65_T, p->H66_T},
                   p->j_dir, data_location, p->Y11==1),
      m_shared_mf(nullptr)
{
    field_amrex::p = p;
    mf = make_mf(p, p->ncomp, &mf, static_cast<int>(data_location));

    BCRecs.resize(p->nlevs);
    for (auto& bc_rec : BCRecs)
        bc_rec.resize(p->ncomp);

    p->register_field(this);
}

field_amrex::~field_amrex()
{
    if (p && m_shared_mf == nullptr && !mf.empty())
        p->deregister_mf(&mf);

    if (p != nullptr)
    p->deregister_field(this);
}

void field_amrex::extend_levels(int new_nlevs)
{
    const int old_nlevs = BCRecs.size();
    if (new_nlevs <= old_nlevs) return;

    // Extend BCRecs with new empty vectors
    BCRecs.resize(new_nlevs);
    for (int lev = old_nlevs; lev < new_nlevs; ++lev)
        BCRecs[lev].resize(p->ncomp);

    // Extend m_alias (view mode only)
    if (m_shared_mf)
    {
        m_alias.resize(new_nlevs);
        for (int lev = old_nlevs; lev < new_nlevs; ++lev)
            m_alias[lev] = amrex::MultiFab((*m_shared_mf)[lev], amrex::make_alias, m_comp, 1);
    }
}

void field_amrex::CopyFrom(const field& src)
{
    const field_amrex* src_amrex = dynamic_cast<const field_amrex*>(&src);
    if (!src_amrex)
        throw std::runtime_error("field_amrex::CopyFrom: source field is not of type field_amrex");

    if (p->nlevs != src_amrex->p->nlevs)
        throw std::runtime_error("field_amrex::CopyFrom: source and destination fields have different number of levels");

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        get_mf(lev).copy(src_amrex->get_mf_const(lev), 0, 0, 1, 0, 0);
    }
}


void field_amrex::rebuild_alias_level(int lev)
{
    // Always invalidate the Array4 cache — the underlying MultiFab storage at lev
    // has been replaced (owning mode) or is about to be re-aliased (view mode).
    m_cached_level         = -1;
    m_cached_mfi_idx       = -1;
    m_cached_const_level   = -1;
    m_cached_const_mfi_idx = -1;

    if (!m_shared_mf) return;  // owning mode: nothing more to do
    if (lev < 0 || lev >= (int)m_alias.size()) return;
    m_alias[lev] = amrex::MultiFab((*m_shared_mf)[lev], amrex::make_alias, m_comp, 1);
}

// ---------------------------------------------------------------------------
// View constructor — non-owning view into shared_mf at component comp
// ---------------------------------------------------------------------------
field_amrex::field_amrex(lexer* p, amrex::Vector<amrex::MultiFab>* shared_mf, int comp,
                         amrex_bc_func::DataLocation data_location)
    : const_params{{p->bcside1, p->bcside4, p->bcside3, p->bcside2, p->bcside5, p->bcside6},
                   {p->H61_T, p->H64_T, p->H63_T, p->H62_T, p->H65_T, p->H66_T},
                   bool(p->j_dir), data_location, p->Y11==1},
      m_shared_mf(shared_mf)
{
    field_amrex::p = p;
    m_comp = comp;
    // mf stays empty — storage is owned by the caller via shared_mf

    BCRecs.resize(p->nlevs);
    for (auto& bc_rec : BCRecs)
        bc_rec.resize(p->ncomp);

    // Build 1-component aliases for GetMultiFab() and get_alias()
    m_alias.resize(p->nlevs);
    for (int lev = 0; lev < m_alias.size(); ++lev)
        m_alias[lev] = amrex::MultiFab((*shared_mf)[lev], amrex::make_alias, comp, 1);

    p->register_field(this);
}

// ---------------------------------------------------------------------------
// setVal
// ---------------------------------------------------------------------------
void field_amrex::setVal(double val, bool includeGhost)
{
    LEVEL_LOOP
    {
        const auto ngrow = includeGhost ? get_mf(p->level).nGrowVect() : amrex::IntVect{0};
        get_mf(p->level).setVal(val, 0, 1, ngrow);
    }
}

// ---------------------------------------------------------------------------
// FillBoundary (MPI ghost exchange, component-restricted)
// ---------------------------------------------------------------------------
void field_amrex::FillBoundary()
{
    LEVEL_LOOP
    {
        get_mf(p->level).FillBoundary(0, 1,
                                       p->amrex_geometry[p->level].periodicity());
    }
}

// ---------------------------------------------------------------------------
// FillBoundaryBatch — one FillBoundary call for a component range
// ---------------------------------------------------------------------------
void field_amrex::FillBoundaryBatch(lexer* p,
                                     amrex::Vector<amrex::MultiFab>& shared_mf,
                                     int scomp, int ncomp)
{
    LEVEL_LOOP
    {
        shared_mf[p->level].FillBoundary(scomp, ncomp,
                                          p->amrex_geometry[p->level].periodicity());
    }
}

// ---------------------------------------------------------------------------
// FillDomainBoundaryValue
// ---------------------------------------------------------------------------
void field_amrex::FillDomainBoundaryValue(double value, int dir, bool high)
{
    LEVEL_LOOP
    {
        amrex::Box dom = p->amrex_geometry[p->level].Domain();
        TILE_LOOP
        {
            // Valid box of the current FAB: index() is by definition its BoxArray subscript.
            const amrex::Box validbox = p->amr_cell_mf[p->level].boxArray()[p->amr_fab_mfi_idx];
            amrex::Box gbx = validbox;
            bool apply = false;

            amrex::IntVect ng(p->margin);
            if (!const_params.y_dimension_exists) ng[1] = 0;
            gbx.grow(ng);

            if ((validbox.smallEnd(0) == dom.smallEnd(0)) && (dir == 0 && !high))
            { gbx.setBig(0, dom.smallEnd(0) - 1); apply = true; }
            if ((validbox.bigEnd(0) == dom.bigEnd(0)) && (dir == 0 && high))
            { gbx.setSmall(0, dom.bigEnd(0) + 1); apply = true; }

            if ((validbox.smallEnd(1) == dom.smallEnd(1)) && (dir == 1 && !high))
            { gbx.setBig(1, dom.smallEnd(1) - 1); apply = true; }
            if ((validbox.bigEnd(1) == dom.bigEnd(1)) && (dir == 1 && high))
            { gbx.setSmall(1, dom.bigEnd(1) + 1); apply = true; }

            if ((validbox.smallEnd(2) == dom.smallEnd(2)) && (dir == 2 && !high))
            { gbx.setBig(2, dom.smallEnd(2) - 1); apply = true; }
            if ((validbox.bigEnd(2) == dom.bigEnd(2)) && (dir == 2 && high))
            { gbx.setSmall(2, dom.bigEnd(2) + 1); apply = true; }

            if (apply)
            {
                auto arr = get_array(p->level, p->amr_local_fab_idx);
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    arr(i, j, k) = value;
                });
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ShiftBigBoundaryFaceInward (static)
// ---------------------------------------------------------------------------
void field_amrex::ShiftBigBoundaryFaceInward(amrex::MultiFab& mf_in,
                                              amrex_bc_func::DataLocation data_location,
                                              const amrex::Geometry& geom)
{
    int dir = -1;
    if      (data_location == amrex_bc_func::DataLocation::FACE_X) dir = 0;
    else if (data_location == amrex_bc_func::DataLocation::FACE_Y) dir = 1;
    else if (data_location == amrex_bc_func::DataLocation::FACE_Z) dir = 2;

    if (dir == -1) return;

    const int domain_hi = geom.Domain().bigEnd(dir);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mf_in); mfi.isValid(); ++mfi)
    {
        const amrex::Box& valid_box = mfi.validbox();

        if (valid_box.bigEnd(dir) != domain_hi) continue;

        const amrex::Box& box = mfi.fabbox();
        amrex::Array4<amrex::Real> const& arr = mf_in.array(mfi);

        int start = domain_hi;
        int end   = box.bigEnd(dir) - 1;

        amrex::Box para_box = box;
        para_box.setSmall(dir, start);
        para_box.setBig(dir, start);

        if (dir == 0)
        {
            amrex::ParallelFor(para_box, mf_in.nComp(),
            [=] AMREX_GPU_DEVICE (int /*i*/, int j, int k, int n)
            {
                for (int i = start; i <= end; ++i)
                    arr(i, j, k, n) = arr(i + 1, j, k, n);
            });
        }
        else if (dir == 1)
        {
            amrex::ParallelFor(para_box, mf_in.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int /*j*/, int k, int n)
            {
                for (int j = start; j <= end; ++j)
                    arr(i, j, k, n) = arr(i, j + 1, k, n);
            });
        }
        else
        {
            amrex::ParallelFor(para_box, mf_in.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/, int n)
            {
                for (int k = start; k <= end; ++k)
                    arr(i, j, k, n) = arr(i, j, k + 1, n);
            });
        }
    }
}

// =========================================================================
// FillCoarseFineNormalGhost — stagger-correct C-F normal-velocity ghost fill.
// FillPatchTwoLevels used cell_cons_interp (cell-centred) on the staggered face
// velocity, which mis-places the C-F normal velocity by ~half a fine cell and
// pulls in the covered coarse face. Here we overwrite the C-F ghosts in the
// field's normal direction with a face-linear interpolation of the coarse face
// velocity (the divergence-preserving prolongation). Domain-boundary ghosts
// (BC functor) and fine-fine ghosts (FillBoundary) are left untouched.
// =========================================================================
void field_amrex::FillCoarseFineNormalGhost()
{
    int dir = -1;
    if      (const_params.data_location == amrex_bc_func::DataLocation::FACE_X) dir = 0;
    else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Y) dir = 1;
    else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Z) dir = 2;
    else return;                       // cell-centred field
    if (p->nlevs <= 1) return;

    const int r = p->ref_vec[dir];
    const amrex::IntVect rv = p->ref_vec;

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        amrex::MultiFab&       fine_mf = GetMultiFab(lev);
        const amrex::Box       dom     = p->amrex_geometry[lev].Domain();
        const amrex::BoxArray& fba     = fine_mf.boxArray();
        const int              ngdir   = fine_mf.nGrow(dir);
        if (ngdir <= 0) continue;

        // coarse velocity on the coarsened-fine layout, with enough ghost in dir
        amrex::BoxArray cfba = amrex::coarsen(fba, rv);
        const int       cng  = ngdir / r + 2;
        amrex::MultiFab cmf(cfba, fine_mf.DistributionMap(), 1, cng);
        cmf.setVal(0.0);
        cmf.ParallelCopy(GetMultiFab(lev-1), 0, 0, 1,
                         amrex::IntVect(0), amrex::IntVect(cng),
                         p->amrex_geometry[lev-1].periodicity());

        for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& vbx = mfi.validbox();
            const amrex::Box& fbx = mfi.fabbox();
            auto       farr = fine_mf.array(mfi);
            const auto carr = cmf.const_array(mfi);

            for (int side = 0; side < 2; ++side)
            {
                amrex::Box gbx = vbx;
                if (side == 0) { gbx.setSmall(dir, vbx.smallEnd(dir)-ngdir); gbx.setBig(dir, vbx.smallEnd(dir)-1); }
                else           { gbx.setSmall(dir, vbx.bigEnd(dir)+1);       gbx.setBig(dir, vbx.bigEnd(dir)+ngdir); }
                gbx &= fbx;
                if (gbx.isEmpty()) continue;

                amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) noexcept
                {
                    const amrex::IntVect iv(i,j,k);
                    if (!dom.contains(iv)) return;   // domain-boundary ghost: BC functor owns it
                    if (fba.contains(iv))  return;   // fine-fine ghost: FillBoundary owns it

                    // face-linear interpolation of the coarse face velocity in `dir`.
                    // coarse face w_coarse(kc) sits at fine-face index r*(kc+1)-1.
                    const int    kf  = iv[dir];
                    const int    kc  = (kf + 1) / r - 1;
                    const double wgt = double((kf + 1) - (kc + 1) * r) / double(r);
                    amrex::IntVect ic  = amrex::coarsen(iv, rv);
                    amrex::IntVect ic2 = ic;
                    ic[dir]  = kc;
                    ic2[dir] = kc + 1;
                    farr(i,j,k) = (1.0 - wgt) * carr(ic) + wgt * carr(ic2);
                });
            }
        }
    }
}

// =========================================================================
// FillCoarseFineCellGhost — matrix-consistent C-F ghost fill for cell-centred
// pressure (gcv 41). FillPatchTwoLevels filled the C-F ring with cell_cons_interp,
// which injects a limited transverse slope built from neighbouring coarse cells. The
// SSAMG matrix row of a fine boundary cell instead couples to exactly ONE coarse cell,
// rescaled to the centre distance d_cf = 0.5*(dx_f+dx_c) (hypre_ssamg_fill.cpp pass 2a).
// Here we overwrite the first fine C-F ghost layer with the two-point normal prolongation
//   p_ghost = p_bnd + (2/(1+r)) * (p_coarse - p_bnd),
// for which (p_ghost - p_bnd)/dx_f == (p_coarse - p_bnd)/d_cf exactly. The corrector's
// fine-spacing gradient then equals the matrix flux, so the assembled operator equals
// G_v^T W G_v across the interface. Domain-boundary ghosts (BC functor) and fine-fine
// ghosts (FillBoundary) are left untouched.
//
// Pass 2 (coarse side): the covered coarse cell C' that an UNcovered coarse cell C reads
// as its C-F neighbour (e.g. the pressure-predictor gradient on the coarse level, which
// has no sole-writer). Prolongate it symmetrically so C's dx_c-spacing gradient reproduces
// the d_cf gradient to the fine average p_favg = average_down(fine) onto C':
//   p_C' = p_C + (2r/(1+r)) * (p_favg - p_C),
// for which (p_C - p_C')/dx_c == (p_C - p_favg)/d_cf. This mirrors the conservative coarse
// C-F coupling (hypre_ssamg_fill.cpp pass 2b). Exact where the C-F interface is axis-
// separated (each covered ring cell has a single uncovered normal neighbour, e.g. a fully
// submerged xy-slab); at a patch corner a covered cell borders uncovered cells on >1 axis
// and only one axis' value can be stored — the last written wins (a cell-ghost limitation;
// full generality would need a sole-writer face correction as cf_velocity_correction does).
// =========================================================================
void field_amrex::FillCoarseFineCellGhost(bool transverse)
{
    // REEF_CF_PROJECTION_GROUP member (3) — the fine-side (2/(1+r)) and coarse-side
    // (2r/(1+r)) prolongations here must match the matrix C-F coupling in
    // hypre_ssamg::amr_cf_coefficients (canonical note there). Change one -> change all.
    // transverse=true (gcv 42, PREDICTOR fields only) additionally reconstructs the coarse
    // value at the fine ghost's TRANSVERSE sub-position instead of using the raw coarse cell,
    // so a field varying transverse to the C-F face keeps its correct normal gradient. Never
    // used for pcorr (would desync from the matrix flux); see FillCoarseFineCellGhost doc.
    if (const_params.data_location != amrex_bc_func::DataLocation::CELL_CENTERED) return;
    if (p->nlevs <= 1) return;

    const amrex::IntVect rv = p->ref_vec;

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        amrex::MultiFab&       fine_mf = GetMultiFab(lev);
        const amrex::Box       dom     = p->amrex_geometry[lev].Domain();
        const amrex::BoxArray& fba     = fine_mf.boxArray();

        // Coarse pressure on the coarsened-fine layout, with a 2-cell ghost so that
        // coarsen(ghost) — which can sit one coarse cell outside the coarsened validbox —
        // is always addressable. ParallelCopy pulls in the (uncovered) coarse-side cells.
        amrex::BoxArray cfba = amrex::coarsen(fba, rv);
        const int       cng  = 2;
        amrex::MultiFab cmf(cfba, fine_mf.DistributionMap(), 1, cng);
        cmf.setVal(0.0);
        cmf.ParallelCopy(GetMultiFab(lev-1), 0, 0, 1,
                         amrex::IntVect(0), amrex::IntVect(cng),
                         p->amrex_geometry[lev-1].periodicity());

        for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& vbx = mfi.validbox();
            const amrex::Box& fbx = mfi.fabbox();
            auto       farr = fine_mf.array(mfi);
            const auto carr = cmf.const_array(mfi);

            for (int dir = 0; dir < 3; ++dir)
            {
                if (dir == 1 && p->j_dir != 1) continue;   // degenerate y (2D x-z run)

                const double         s = 2.0 / (1.0 + double(rv[dir]));
                const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

                for (int side = 0; side < 2; ++side)
                {
                    // First ghost layer outboard of the valid box in `dir`. Transverse
                    // extent = vbx (only face-adjacent ghosts feed the axis-aligned corrector).
                    amrex::Box gbx = vbx;
                    if (side == 0) { gbx.setSmall(dir, vbx.smallEnd(dir)-1); gbx.setBig(dir, vbx.smallEnd(dir)-1); }
                    else           { gbx.setSmall(dir, vbx.bigEnd(dir)+1);   gbx.setBig(dir, vbx.bigEnd(dir)+1);   }
                    gbx &= fbx;
                    if (gbx.isEmpty()) continue;

                    amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) noexcept
                    {
                        const amrex::IntVect g(i,j,k);
                        if (!dom.contains(g)) return;   // domain-boundary ghost: BC functor owns it
                        if (fba.contains(g))  return;   // fine-fine ghost: FillBoundary owns it

                        const amrex::IntVect b  = (side == 0) ? g + e : g - e;   // interior boundary cell
                        const double         pb = farr(b);
                        const amrex::IntVect cc = amrex::coarsen(g, rv);
                        double               pc = carr(cc);

                        // Transverse-linear reconstruction: evaluate the coarse field at the fine
                        // ghost's sub-position within its coarse cell, so a field varying transverse
                        // to `dir` keeps the correct normal gradient (raw carr(cc) is transverse-
                        // piecewise-constant -> spurious normal gradient = the geyser). Central slope
                        // per transverse axis (cmf carries a 2-cell ghost, so cc +/- et is addressable).
                        if (transverse)
                            for (int t = 0; t < 3; ++t)
                            {
                                if (t == dir) continue;
                                if (t == 1 && p->j_dir != 1) continue;
                                const amrex::IntVect et = amrex::IntVect::TheDimensionVector(t);
                                const double slope = 0.5 * (carr(cc + et) - carr(cc - et));
                                const int    sub   = g[t] - cc[t] * rv[t];          // 0 .. r-1
                                const double off   = (double(sub) + 0.5) / double(rv[t]) - 0.5; // (-1/2,1/2)
                                pc += slope * off;
                            }

                        farr(i,j,k) = pb + s * (pc - pb);
                    });
                }
            }
        }
    }

    // --- Pass 2: coarse-side covered ring -------------------------------------------
    // Only the matrix-consistent (pcorr) path needs the coarse-side d_cf prolongation of the
    // covered ring. The transverse predictor path (gcv 42) must NOT touch covered coarse cells:
    // they are owned by average_down / REEF_COVERED_PRESS_RECON (or skipped via skip_covered),
    // and the d_cf overwrite here re-injects a spurious offset (breaks the hydrostatic balance).
    if (transverse) return;

    for (int lev = 0; lev < p->nlevs - 1; ++lev)
    {
        amrex::MultiFab& coarse_mf = GetMultiFab(lev);

        // Fine average of level lev+1 onto this coarse layout. average_down writes only the
        // covered cells (uncovered coarse cells have no fine coverage), so p_favg is defined
        // exactly where we need it (the covered ring).
        amrex::MultiFab favg(p->amrex_box_array[lev], p->amrex_distribution_mapping[lev], 1, 0);
        favg.setVal(0.0);
        amrex::average_down(GetMultiFab(lev+1), favg, 0, 1, rv);

        // Covered mask with a 1-cell ghost so the normal neighbour's covered/uncovered state
        // is known across box boundaries (amr_cell_mf: 1 uncovered, 0 covered). Domain-exterior
        // ghosts default to 1 (uncovered) — covered cells never touch the domain boundary.
        amrex::iMultiFab maskg(p->amr_cell_mf[lev].boxArray(),
                               p->amr_cell_mf[lev].DistributionMap(), 1, 1);
        maskg.setVal(1);
        amrex::iMultiFab::Copy(maskg, p->amr_cell_mf[lev], 0, 0, 1, 0);
        maskg.FillBoundary(p->amrex_geometry[lev].periodicity());

        for (amrex::MFIter mfi(coarse_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& vbx = mfi.validbox();
            auto       parr = coarse_mf.array(mfi);
            const auto favgr = favg.const_array(mfi);
            const auto maskr = maskg.const_array(mfi);

            amrex::LoopOnCpu(vbx, [&] (int i, int j, int k) noexcept
            {
                if (maskr(i,j,k) != 0) return;   // uncovered: not a covered ring cell

                const amrex::IntVect c(i,j,k);
                for (int dir = 0; dir < 3; ++dir)
                {
                    if (dir == 1 && p->j_dir != 1) continue;

                    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
                    const double s_c = 2.0 * double(rv[dir]) / (1.0 + double(rv[dir]));

                    // Uncovered normal neighbour on either side across the C-F interface.
                    const bool up = (maskr(c + e) != 0);
                    const bool dn = (maskr(c - e) != 0);
                    if (!up && !dn) continue;

                    const double pC   = up ? parr(c + e) : parr(c - e);
                    const double pfav = favgr(c);
                    parr(i,j,k) = pC + s_c * (pfav - pC);
                }
            });
        }
    }
}

// =========================================================================
// average_down_level — stagger-correct fine->coarse average down (lev+1 -> lev).
// Cell-centred fields (data_location == CELL_CENTERED) use the plain cell-centred
// amrex::average_down. FACE_X/Y/Z fields (data_location 1/2/3) hold the staggered
// velocity in a cell-centred MultiFab with vel(cell i) == the high face == face(i+e),
// so a cell-centred average_down would mix face values across the C-F interface and
// inject divergence. Instead we convert to a face-typed temporary, reflux with
// amrex::average_down_faces, and copy the synced coarse faces back. Only the field's
// own normal direction is processed.
// =========================================================================
void field_amrex::average_down_level(lexer* p, int lev)
{
    if (lev < 0 || lev + 1 >= p->nlevs) return;

    const int flev = lev + 1;
    const int clev = lev;

    if (const_params.data_location == amrex_bc_func::DataLocation::CELL_CENTERED)
    {
        amrex::average_down(GetMultiFab(flev), GetMultiFab(clev), 0, 1, p->ref_vec);
        return;
    }

    int dir = -1;
    if      (const_params.data_location == amrex_bc_func::DataLocation::FACE_X) dir = 0;
    else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Y) dir = 1;
    else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Z) dir = 2;
    else return;

    if (dir == 1 && p->j_dir != 1) return;   // degenerate y (2D x-z run): FACE_Y unused

    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

    amrex::MultiFab ffine(amrex::convert(p->amrex_box_array[flev], e),
                          p->amrex_distribution_mapping[flev], 1, 0);
    amrex::MultiFab fcrse(amrex::convert(p->amrex_box_array[clev], e),
                          p->amrex_distribution_mapping[clev], 1, 0);

    // staggered cell-centred velocity -> face MultiFab (face f = vel(f-e))
    auto fill_face = [&] (amrex::MultiFab& fmf, int amrlev)
    {
        auto& v_mf = GetMultiFab(amrlev);
        for (amrex::MFIter mfi(fmf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& fbx = mfi.validbox();
            auto       f  = fmf.array(mfi);
            const auto vv = v_mf.const_array(mfi);
            amrex::LoopOnCpu(fbx, [&] (int i, int j, int k) noexcept
            { f(i,j,k) = vv(i-e[0], j-e[1], k-e[2]); });
        }
    };
    fill_face(ffine, flev);
    fill_face(fcrse, clev);   // pre-fill so uncovered coarse faces survive

    amrex::average_down_faces(ffine, fcrse, p->ref_vec, p->amrex_geometry[clev]);

    // synced coarse faces -> cell-centred storage (vel(i) = face(i+e))
    auto& v_mf = GetMultiFab(clev);
    for (amrex::MFIter mfi(v_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        auto       vv = v_mf.array(mfi);
        const auto f  = fcrse.const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        { vv(i,j,k) = f(i+e[0], j+e[1], k+e[2]); });
    }
}

#endif