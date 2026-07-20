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

#include "ArrayWrapper3D.h"
#include "ArrayWrapper3D_imp.h"
#include "lexer.h"
#include <algorithm>

ArrayWrapper3D::ArrayWrapper3D(lexer *pp, DataLocation _data_location) : p(pp), data_location(_data_location)
{
}

#if USE_AMREX
ArrayWrapper3D::ArrayWrapper3D(lexer *pp, amrex::Vector<amrex::iMultiFab> *shared, int comp)
    : p(pp), m_shared(shared), m_comp(comp)
{
}
#endif

ArrayWrapper3D::~ArrayWrapper3D()
{
    #if USE_AMREX
    if (p && !m_shared && !data.empty())
        p->deregister_imf(&data);
    #endif
}

int ArrayWrapper3D::kz() const noexcept
{
    return data_location == DataLocation::NODE_Z ? p->kmaxF : p->kmax;
}


void ArrayWrapper3D::resize(int default_value)
{
    #if USE_AMREX
    if (m_shared) return; // view mode: shared MultiFab is resized by the owner
    if (data.empty())
        p->register_imf(&data, 1);
    // NODE_Z needs a z-nodal BoxArray: one z-plane more
    // than there are cells. amrex::convert shares the underlying box list, so
    // box count, ordering and the DistributionMap are preserved and the
    // existing MFIter walks it unchanged.
    //
    // register_imf carries no data location yet, so the registry would redefine
    // this on the cell-centred BoxArray at the first regrid. Nothing addresses
    // location 7 under AMReX until that registration becomes location-aware.
    data.resize(p->nlevs);
    LEVEL_LOOP
    {
        const amrex::BoxArray ba = (data_location == DataLocation::NODE_Z)
            ? amrex::convert(p->amrex_box_array[p->level], amrex::IntVect(AMREX_D_DECL(0,0,1)))
            : p->amrex_box_array[p->level];

        data[p->level].define(ba, p->amrex_distribution_mapping[p->level], 1, p->margin);
        data[p->level].setVal(default_value);
    }
    #else
    // Single level: one flat array, sized from the same lexer metrics the IJK
    // (or, for NODE_Z, FIJK) macro reads. grid::assign_margin sets
    // imax/jmax/kmax/kmaxF and imin/jmin/kmin together, so all are final by the
    // time resize runs.
    //
    // The NODE_Z slack plane matches the imax*jmax*(kmax+2) allocation this
    // replaced in driver_makegrid_sigma.cpp: the stride is kmaxF = kmax+1, and
    // the forward-stencil macros (FIJKp3/p4) reach past the last in-stride slot
    // in the final column. See field7 for the same reasoning.
    const std::size_t plane = static_cast<std::size_t>(p->imax)*p->jmax;
    const std::size_t slack = (data_location == DataLocation::NODE_Z) ? plane : 0;
    data.resize(plane*kz() + slack, default_value);
    cache_addressing();
    #endif
}

void ArrayWrapper3D::setVal(int val, bool includeGhost)
{
    #if USE_AMREX
    LEVEL_LOOP
    {
        const int ng = includeGhost ? p->margin : 0;
        GetMultiFab(p->level).setVal(val, m_comp, 1, ng);
    }
    #else
    if(includeGhost)
    {
        std::fill(data.begin(), data.end(), val);
    }
    else if(data_location == DataLocation::NODE_Z)
    {
        // FBASELOOP, not LOOP: LOOP stops at KMAX_LOOP and would leave the top
        // node plane untouched, and its PCHECK reads the IJK-strided flag4.
        int i,j,k;
        FBASELOOP
        {
            operator()(i,j,k) = val;
        }
    }
    else
    {
        int i,j,k;
        LOOP
        {
            operator()(i,j,k) = val;
        }
    }
    #endif
}

ArrayWrapper3D::operator int *()
{
    #if USE_AMREX
    return GetMultiFab(p->level).atLocalIdx(p->amr_local_fab_idx).dataPtr(m_comp);
    #else
    return data.data(); // unshifted: callers index this with IJK, not with i/j/k
    #endif
}

#if USE_AMREX
void ArrayWrapper3D::fillBoundary()
{
    LEVEL_LOOP
    {
        // NODE_Z valid regions overlap on z-split box seams
        // — the shared plane is valid in both neighbours and FillBoundary only
        // fills ghosts, it does not arbitrate between two valid copies.
        // OverrideSync picks the canonical owner first; this is what
        // gcx_parax7co does for the legacy flat arrays.
        if (data_location == DataLocation::NODE_Z)
            GetMultiFab(p->level).OverrideSync(m_comp, 1, p->amrex_geometry[p->level].periodicity());

        GetMultiFab(p->level).FillBoundary(m_comp, 1);
    }
}

void ArrayWrapper3D::FillBoundaryBatch(lexer* p, amrex::Vector<amrex::iMultiFab>& shared,
                                         int scomp, int ncomp)
{
    for (int lev = 0; lev < p->nlevs; ++lev)
        shared[lev].FillBoundary(scomp, ncomp);
}

void ArrayWrapper3D::fillHigherLevels()
{
    const amrex::IntVect ratio = p->ref_vec;
    const int ratio_x = ratio[0];
    const int ratio_y = ratio[1];
    const int ratio_z = ratio[2];

    int dir = -1;
    if (data_location == DataLocation::FACE_X) dir = 0;
    else if (data_location == DataLocation::FACE_Y) dir = 1;
    else if (data_location == DataLocation::FACE_Z) dir = 2;

    for (int lev = 1; lev < p->nlevs; ++lev)
    {
        const auto& crse_mf = GetMultiFab(lev-1);
        auto& fine_mf = GetMultiFab(lev);

        // Coarse data on the fine layout, carrying enough ghost cells that the COARSENED fine
        // ghost ring is covered. The fine ghost ring at a coarse-fine interface has no
        // same-level neighbour, so the FillBoundary() below cannot fill it; injecting the coarse
        // value here is what gives those C-F ghosts the proper (water) flag instead of the
        // resize(-1) default. Without it poisson_pcorr reads flag<0 across the interface, zeros
        // the C-F face, and destroys the amr_cf coupling (-> projection divergence blow-up).
        const int ng_f = fine_mf.nGrow();
        int rmax = ratio_x;
        if (ratio_y > rmax) rmax = ratio_y;
        if (ratio_z > rmax) rmax = ratio_z;
        const int ng_c = (ng_f + rmax - 1) / rmax;

        amrex::BoxArray coarsened_fine_ba = amrex::coarsen(fine_mf.boxArray(), ratio);

        amrex::iMultiFab coarse_on_fine_layout(coarsened_fine_ba, fine_mf.DistributionMap(), 1, ng_c);
        coarse_on_fine_layout.ParallelCopy(crse_mf, m_comp, 0, 1,
                                           amrex::IntVect(0), amrex::IntVect(ng_c),
                                           p->amrex_geometry[lev-1].periodicity());

        // Only fill cells inside the physical domain: valid cells and C-F interface ghosts.
        // Domain-exterior ghosts (true walls) must keep their OBJ_FLAG.
        const amrex::Box fine_domain = p->amrex_geometry[lev].Domain();

        for (amrex::MFIter mfi(fine_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& fine_valid_box = mfi.validbox();
            const amrex::Box  fill_box = amrex::grow(fine_valid_box, ng_f);
            auto const& fine_arr = fine_mf.array(mfi);
            auto const& crse_arr = coarse_on_fine_layout.const_array(mfi);

            amrex::ParallelFor(fill_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (!fine_domain.contains(amrex::IntVect(i, j, k))) return;
                const int ic = amrex::coarsen(i, ratio_x);
                const int jc = amrex::coarsen(j, ratio_y);
                const int kc = amrex::coarsen(k, ratio_z);
                fine_arr(i, j, k, m_comp) = crse_arr(ic, jc, kc, 0);
            });

            if(dir != -1)
            {
                if(fine_valid_box.bigEnd(dir) == p->amrex_geometry[lev].Domain().bigEnd(dir))
                {
                    amrex::Box end_box = fine_valid_box;
                    end_box.setBig(dir, p->amrex_geometry[lev].Domain().bigEnd(dir)-1);
                    end_box.setSmall(dir, p->amrex_geometry[lev].Domain().bigEnd(dir)-2);

                    int ii = 0, jj = 0, kk = 0;
                    if(dir==0) ii = -1;
                    else if(dir==1) jj = -1;
                    else if(dir==2) kk = -1;
                    amrex::ParallelFor(end_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        fine_arr(i, j, k, m_comp) = fine_arr(i+ii, j+jj, k+kk, m_comp);
                    });
                }
            }
        }

        fine_mf.FillBoundary();
    }
}

#else

/*!
    * @brief Precomputes the flat addressing used by operator() and operator[].
    *
    * The origin and both strides are folded into m_base, so that
    *
    *   data[(i-imin)*jmax*kmax + (j-jmin)*kmax + (k-kmin)]  ==  m_base[i*m_js + j*m_ks + k]
    *
    * and operator() needs no lexer dereference at all. That is what lets the
    * index be shared with neighbouring accesses instead of being rebuilt from
    * p on every one — the wrapper's own p is a different pointer to the
    * compiler than the p at the call site, so nothing computed through it can
    * be reused.
    *
    * m_base is not a pointer formed before the array: imin/jmin/kmin are all
    * -margin, so it lands at data.data() + margin*(m_js + m_ks + 1).
    *
    * The strides are long deliberately. An int store through the reference
    * returned by operator() may alias an int member under TBAA, which forces a
    * reload of the strides on every iteration of a writing loop; as longs they
    * stay in registers and the address strength-reduces to an induction
    * variable.
*/
void ArrayWrapper3D::cache_addressing() noexcept
{
    m_ks   = kz();
    m_js   = static_cast<long>(p->jmax) * kz();
    m_base = data.data() - p->imin*m_js - p->jmin*m_ks - p->kmin;
}
#endif
