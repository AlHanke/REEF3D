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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include "hypre_ssamg.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "field.h"
#include "fieldint4.h"
#include <HYPRE_parcsr_mv.h>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <set>
#include <vector>

// Two cheap consistency pre-checks on the assembled pressure-correction operator A.
// Both are matrix-only: they apply A (via HYPRE's own matvec, so the C-F non-stencil
// couplings and the solid identity rows are exercised exactly as solved) to test
// vectors and inspect the result. Nothing here touches a->M directly, so the checks
// see precisely the operator handed to the solver.
//
//   1. rowsum   A*1 : 0 on a clean interior fluid row (the stencil annihilates a
//                     constant); on a boundary row it equals the face transmissibilities
//                     that poisson_pcorr zeroed into the RHS but left on the diagonal.
//                     A non-zero interior value flags a diagonal/off-diagonal imbalance
//                     or a miscounted C-F coupling.
//   2. symmetry      global  yT A x - xT A y  (0 for a symmetric operator). Away from the
//                     C-F interface the operator must be symmetric; the conservative C-F
//                     coupling is deliberately not. On the ParCSR (nlevs>1) path the
//                     per-row  Ax - A^T x  localises the asymmetry: it must vanish on
//                     interior rows and be confined to the C-F ring.

void hypre_ssamg::validate_operator(lexer* p, fdm* a, ghostcell* pgc)
{
    const int rank = p->mpirank;

    // --- shared helpers --------------------------------------------------------
    auto vcreate = [&](HYPRE_SStructVector& v)
    {
        HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &v);
        HYPRE_SStructVectorSetObjectType(v, object_type);
        HYPRE_SStructVectorInitialize(v);
    };

    // y = A x  (transpose: y = A^T x, ParCSR path only). For the ParCSR object the
    // structured side of y is stale after the par matvec, so gather it back -- the same
    // refresh solve() performs after the par solve.
    auto matvec = [&](HYPRE_SStructVector xv, HYPRE_SStructVector yv, bool transpose)
    {
        if (p->nlevs > 1)
        {
            HYPRE_ParCSRMatrix pA; HYPRE_ParVector px, py;
            HYPRE_SStructMatrixGetObject(A,  (void**) &pA);
            HYPRE_SStructVectorGetObject(xv, (void**) &px);
            HYPRE_SStructVectorGetObject(yv, (void**) &py);
            if (transpose) HYPRE_ParCSRMatrixMatvecT(1.0, pA, px, 0.0, py);
            else           HYPRE_ParCSRMatrixMatvec (1.0, pA, px, 0.0, py);
            HYPRE_SStructVectorGather(yv);
        }
        else
        {
            HYPRE_SStructMatrixMatvec(1.0, A, xv, 0.0, yv);
        }
    };

    auto report = [&](const char* tag, double interior, double boundary, double cf)
    {
        const double gi = pgc->globalmax(interior);
        const double gb = pgc->globalmax(boundary);
        const double gc = pgc->globalmax(cf);
        if (rank == 0)
            std::cout << "  [matcheck] " << tag
                      << "  interior=" << gi
                      << "  boundary=" << gb
                      << "  C-F=" << gc << std::endl;
    };

    std::mt19937 rng(12345u + unsigned(rank));
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

#if USE_AMREX
    // C-F "from" cells: every cell whose stencil row was edited by amr_cf_coefficients.
    // Both interface directions are recorded, so the from-cells cover both sides of a face.
    std::set<std::array<int,4>> cfset;
    for (const auto& L : cf_links)
        cfset.insert({L.from_part, L.from_ijk[0], L.from_ijk[1], L.from_ijk[2]});

    auto is_cf = [&](int lev, int i, int j, int k)
    { return cfset.count({lev, i, j, k}) != 0; };

    // A row is a boundary row if a stencil face was dropped: either a neighbour is
    // non-fluid (flag4<0: air/inflow/outflow/solid) or the cell sits on a domain edge.
    // The domain-edge test matters because lateral-wall ghost flag4 can be >=0 (mirror),
    // so the flag check alone would mis-bucket wall cells as interior.
    auto is_boundary = [&](int lev, const amrex::Array4<const int>& fa,
                           int i, int j, int k)
    {
        const amrex::Box dom = p->amrex_geometry[lev].Domain();
        if (i == dom.smallEnd(0) || i == dom.bigEnd(0)) return true;
        if (k == dom.smallEnd(2) || k == dom.bigEnd(2)) return true;
        if (p->j_dir && (j == dom.smallEnd(1) || j == dom.bigEnd(1))) return true;
        if (fa(i-1,j,k) < 0 || fa(i+1,j,k) < 0
         || fa(i,j,k-1) < 0 || fa(i,j,k+1) < 0
         || (p->j_dir && (fa(i,j-1,k) < 0 || fa(i,j+1,k) < 0))) return true;
        return false;
    };

    // Set v on every fluid cell from gen() (sequential draws / constant), 0 elsewhere.
    auto fill = [&](HYPRE_SStructVector v, const std::function<double()>& gen)
    {
        for (int lev = 0; lev < p->nlevs; ++lev)
        {
            const auto& flag4_mf = p->flag4.GetMultiFab(lev);
            for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
                int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
                const auto fa = flag4_mf.const_array(mfi);

                std::vector<double> vals(bx.numPts());
                int cnt = 0;
                for (int kk = lo[2]; kk <= hi[2]; ++kk)
                for (int jj = lo[1]; jj <= hi[1]; ++jj)
                for (int ii = lo[0]; ii <= hi[0]; ++ii)
                    vals[cnt++] = (fa(ii,jj,kk) >= AIR_FLAG) ? gen() : 0.0;

                HYPRE_SStructVectorSetBoxValues(v, lev, lo, hi, variable, vals.data());
            }
        }
        HYPRE_SStructVectorAssemble(v);
    };

    // ---- check 1: rowsum  A*1 -------------------------------------------------
    {
        HYPRE_SStructVector xv, yv;
        vcreate(xv); vcreate(yv);
        HYPRE_SStructVectorAssemble(yv);          // zero output, ready for GetObject
        fill(xv, []{ return 1.0; });
        matvec(xv, yv, false);

        double r_int = 0.0, r_bnd = 0.0, r_cf = 0.0;
        // worst C-F row (rowsum must be ~0 for a conservative C-F coupling); flag whether the
        // worst cell is ALSO a wall/boundary cell (a wall-corner), which is where the Neumann
        // wall fold and the amr_cf diagonal edit both touch M.p and may stop netting to zero.
        double wcf = 0.0; int wcf_lev = -1, wcf_ijk[3] = {-1,-1,-1}, wcf_flag = 0; bool wcf_corner = false;
        for (int lev = 0; lev < p->nlevs; ++lev)
        {
            const auto& flag4_mf = p->flag4.GetMultiFab(lev);
            for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
                int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
                const auto fa = flag4_mf.const_array(mfi);

                std::vector<double> y(bx.numPts());
                HYPRE_SStructVectorGetBoxValues(yv, lev, lo, hi, variable, y.data());

                int cnt = 0;
                for (int kk = lo[2]; kk <= hi[2]; ++kk)
                for (int jj = lo[1]; jj <= hi[1]; ++jj)
                for (int ii = lo[0]; ii <= hi[0]; ++ii)
                {
                    const double val = std::fabs(y[cnt++]);
                    if (fa(ii,jj,kk) < AIR_FLAG) continue;   // identity row, skip

                    if (is_cf(lev, ii, jj, kk))
                    {
                        r_cf = std::max(r_cf, val);
                        if (val > wcf)
                        {
                            wcf = val; wcf_lev = lev;
                            wcf_ijk[0]=ii; wcf_ijk[1]=jj; wcf_ijk[2]=kk; wcf_flag = fa(ii,jj,kk);
                            wcf_corner = is_boundary(lev, fa, ii, jj, kk);
                        }
                    }
                    else if (is_boundary(lev,fa,ii,jj,kk)) r_bnd = std::max(r_bnd, val);
                    else                                   r_int = std::max(r_int, val);
                }
            }
        }
        report("rowsum  A*1     max|.|", r_int, r_bnd, r_cf);
        const double g_cf = pgc->globalmax(r_cf);
        if (wcf_lev >= 0 && wcf == g_cf)
            std::cout << "  [matcheck] worst C-F rowsum=" << wcf << " at lev=" << wcf_lev
                      << " (" << wcf_ijk[0] << "," << wcf_ijk[1] << "," << wcf_ijk[2] << ")"
                      << "  flag4=" << wcf_flag << "  wall-corner=" << wcf_corner << std::endl;

        HYPRE_SStructVectorDestroy(xv);
        HYPRE_SStructVectorDestroy(yv);
    }

    // ---- check 2: symmetry ----------------------------------------------------
    {
        HYPRE_SStructVector xv, yv, Axv, Ayv;
        vcreate(xv); vcreate(yv); vcreate(Axv); vcreate(Ayv);
        HYPRE_SStructVectorAssemble(Axv);
        HYPRE_SStructVectorAssemble(Ayv);
        fill(xv, [&]{ return uni(rng); });
        fill(yv, [&]{ return uni(rng); });
        matvec(xv, Axv, false);
        matvec(yv, Ayv, false);

        // global bilinear asymmetry  yT(Ax) - xT(Ay)
        double s1 = 0.0, s2 = 0.0;
        for (int lev = 0; lev < p->nlevs; ++lev)
        {
            const auto& flag4_mf = p->flag4.GetMultiFab(lev);
            for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
                int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
                const auto fa = flag4_mf.const_array(mfi);
                const int nc = bx.numPts();

                std::vector<double> xx(nc), yy(nc), Ax(nc), Ay(nc);
                HYPRE_SStructVectorGetBoxValues(xv,  lev, lo, hi, variable, xx.data());
                HYPRE_SStructVectorGetBoxValues(yv,  lev, lo, hi, variable, yy.data());
                HYPRE_SStructVectorGetBoxValues(Axv, lev, lo, hi, variable, Ax.data());
                HYPRE_SStructVectorGetBoxValues(Ayv, lev, lo, hi, variable, Ay.data());

                int cnt = 0;
                for (int kk = lo[2]; kk <= hi[2]; ++kk)
                for (int jj = lo[1]; jj <= hi[1]; ++jj)
                for (int ii = lo[0]; ii <= hi[0]; ++ii)
                {
                    if (fa(ii,jj,kk) >= AIR_FLAG)
                    {
                        s1 += yy[cnt] * Ax[cnt];
                        s2 += xx[cnt] * Ay[cnt];
                    }
                    ++cnt;
                }
            }
        }
        const double gs1 = pgc->globalsum(s1);
        const double gs2 = pgc->globalsum(s2);
        const double denom = std::max(std::fabs(gs1), std::max(std::fabs(gs2), 1e-300));
        if (rank == 0)
            std::cout << "  [matcheck] symmetry  yT A x - xT A y = " << (gs1 - gs2)
                      << "  (rel " << std::fabs(gs1 - gs2) / denom << ")" << std::endl;

        // per-row localisation  Ax - A^T x  (ParCSR transpose only)
        if (p->nlevs > 1)
        {
            HYPRE_SStructVector ATx;
            vcreate(ATx);
            HYPRE_SStructVectorAssemble(ATx);
            matvec(xv, ATx, true);

            double d_int = 0.0, d_bnd = 0.0, d_cf = 0.0;
            for (int lev = 0; lev < p->nlevs; ++lev)
            {
                const auto& flag4_mf = p->flag4.GetMultiFab(lev);
                for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.validbox();
                    int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
                    int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
                    const auto fa = flag4_mf.const_array(mfi);
                    const int nc = bx.numPts();

                    std::vector<double> Ax(nc), ATxv(nc);
                    HYPRE_SStructVectorGetBoxValues(Axv, lev, lo, hi, variable, Ax.data());
                    HYPRE_SStructVectorGetBoxValues(ATx, lev, lo, hi, variable, ATxv.data());

                    int cnt = 0;
                    for (int kk = lo[2]; kk <= hi[2]; ++kk)
                    for (int jj = lo[1]; jj <= hi[1]; ++jj)
                    for (int ii = lo[0]; ii <= hi[0]; ++ii)
                    {
                        const double d = std::fabs(Ax[cnt] - ATxv[cnt]);
                        ++cnt;
                        if (fa(ii,jj,kk) < AIR_FLAG) continue;

                        if (is_cf(lev, ii, jj, kk))            d_cf  = std::max(d_cf,  d);
                        else if (is_boundary(lev,fa,ii,jj,kk)) d_bnd = std::max(d_bnd, d);
                        else                                   d_int = std::max(d_int, d);
                    }
                }
            }
            report("symmetry  Ax-A^T x max|.|", d_int, d_bnd, d_cf);

            HYPRE_SStructVectorDestroy(ATx);
        }

        HYPRE_SStructVectorDestroy(xv);
        HYPRE_SStructVectorDestroy(yv);
        HYPRE_SStructVectorDestroy(Axv);
        HYPRE_SStructVectorDestroy(Ayv);
    }

#else
    // Single-part (no AMReX): one box [ilower,iupper], no coarse-fine couplings, so the
    // operator is symmetric and the rowsum splits into interior vs boundary only.
    const int nx = p->knox, ny = p->knoy, nz = p->knoz;
    const int nc = nx * ny * nz;

    auto fill = [&](HYPRE_SStructVector v, const std::function<double()>& gen)
    {
        std::vector<double> vals(nc);
        count = 0;
        KJILOOP
        {
            PFLUIDCHECK vals[count] = gen();
            SFLUIDCHECK vals[count] = 0.0;
            ++count;
        }
        HYPRE_SStructVectorSetBoxValues(v, 0, ilower, iupper, variable, vals.data());
        HYPRE_SStructVectorAssemble(v);
    };

    // ---- check 1: rowsum  A*1 -------------------------------------------------
    {
        HYPRE_SStructVector xv, yv;
        vcreate(xv); vcreate(yv);
        HYPRE_SStructVectorAssemble(yv);
        fill(xv, []{ return 1.0; });
        matvec(xv, yv, false);

        std::vector<double> y(nc);
        HYPRE_SStructVectorGetBoxValues(yv, 0, ilower, iupper, variable, y.data());

        double r_int = 0.0, r_bnd = 0.0;
        count = 0;
        KJILOOP
        {
            PFLUIDCHECK
            {
                const double val = std::fabs(y[count]);
                const bool boundary =
                       p->flag4(i-1,j,k) < 0 || p->flag4(i+1,j,k) < 0
                    || p->flag4(i,j,k-1) < 0 || p->flag4(i,j,k+1) < 0
                    || (p->j_dir && (p->flag4(i,j-1,k) < 0 || p->flag4(i,j+1,k) < 0));
                if (boundary) r_bnd = std::max(r_bnd, val);
                else          r_int = std::max(r_int, val);
            }
            ++count;
        }
        report("rowsum  A*1     max|.|", r_int, r_bnd, 0.0);

        HYPRE_SStructVectorDestroy(xv);
        HYPRE_SStructVectorDestroy(yv);
    }

    // ---- check 2: symmetry (global bilinear) ----------------------------------
    {
        HYPRE_SStructVector xv, yv, Axv, Ayv;
        vcreate(xv); vcreate(yv); vcreate(Axv); vcreate(Ayv);
        HYPRE_SStructVectorAssemble(Axv);
        HYPRE_SStructVectorAssemble(Ayv);
        fill(xv, [&]{ return uni(rng); });
        fill(yv, [&]{ return uni(rng); });
        matvec(xv, Axv, false);
        matvec(yv, Ayv, false);

        std::vector<double> xx(nc), yy(nc), Ax(nc), Ay(nc);
        HYPRE_SStructVectorGetBoxValues(xv,  0, ilower, iupper, variable, xx.data());
        HYPRE_SStructVectorGetBoxValues(yv,  0, ilower, iupper, variable, yy.data());
        HYPRE_SStructVectorGetBoxValues(Axv, 0, ilower, iupper, variable, Ax.data());
        HYPRE_SStructVectorGetBoxValues(Ayv, 0, ilower, iupper, variable, Ay.data());

        double s1 = 0.0, s2 = 0.0;
        count = 0;
        KJILOOP
        {
            PFLUIDCHECK { s1 += yy[count] * Ax[count]; s2 += xx[count] * Ay[count]; }
            ++count;
        }
        const double gs1 = pgc->globalsum(s1);
        const double gs2 = pgc->globalsum(s2);
        const double denom = std::max(std::fabs(gs1), std::max(std::fabs(gs2), 1e-300));
        if (rank == 0)
            std::cout << "  [matcheck] symmetry  yT A x - xT A y = " << (gs1 - gs2)
                      << "  (rel " << std::fabs(gs1 - gs2) / denom << ")" << std::endl;

        HYPRE_SStructVectorDestroy(xv);
        HYPRE_SStructVectorDestroy(yv);
        HYPRE_SStructVectorDestroy(Axv);
        HYPRE_SStructVectorDestroy(Ayv);
    }
#endif
}

// out = A * in on the assembled operator. Fluid rows (flag4>=AIR_FLAG) carry the field
// value; non-fluid rows are zeroed (identity rows then return 0). The matrix A still holds
// the last fill_matrix4 assembly (it is created/destroyed with the object, not per solve).
void hypre_ssamg::matvec_into(lexer* p, fdm* a, ghostcell* pgc, field& out, field& in)
{
    HYPRE_SStructVector xv, yv;
    auto vcreate = [&](HYPRE_SStructVector& v)
    {
        HYPRE_SStructVectorCreate(pgc->mpi_comm, grid, &v);
        HYPRE_SStructVectorSetObjectType(v, object_type);
        HYPRE_SStructVectorInitialize(v);
    };
    vcreate(xv); vcreate(yv);

#if USE_AMREX
    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        const auto& flag4_mf = p->flag4.GetMultiFab(lev);
        const auto& in_mf    = in.GetMultiFab(lev);
        for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};
            const auto fa = flag4_mf.const_array(mfi);
            const auto ia = in_mf.const_array(mfi);

            std::vector<double> vals(bx.numPts());
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
                vals[cnt++] = (fa(ii,jj,kk) >= AIR_FLAG) ? ia(ii,jj,kk) : 0.0;

            HYPRE_SStructVectorSetBoxValues(xv, lev, lo, hi, variable, vals.data());
        }
    }
    HYPRE_SStructVectorAssemble(xv);
    HYPRE_SStructVectorAssemble(yv);

    if (p->nlevs > 1)
    {
        HYPRE_ParCSRMatrix pA; HYPRE_ParVector px, py;
        HYPRE_SStructMatrixGetObject(A,  (void**) &pA);
        HYPRE_SStructVectorGetObject(xv, (void**) &px);
        HYPRE_SStructVectorGetObject(yv, (void**) &py);
        HYPRE_ParCSRMatrixMatvec(1.0, pA, px, 0.0, py);
        HYPRE_SStructVectorGather(yv);
    }
    else
        HYPRE_SStructMatrixMatvec(1.0, A, xv, 0.0, yv);

    for (int lev = 0; lev < p->nlevs; ++lev)
    {
        const auto& flag4_mf = p->flag4.GetMultiFab(lev);
        auto&       out_mf   = out.GetMultiFab(lev);
        for (amrex::MFIter mfi(flag4_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            int lo[3] = {bx.smallEnd(0), bx.smallEnd(1), bx.smallEnd(2)};
            int hi[3] = {bx.bigEnd(0),   bx.bigEnd(1),   bx.bigEnd(2)};

            std::vector<double> vals(bx.numPts());
            HYPRE_SStructVectorGetBoxValues(yv, lev, lo, hi, variable, vals.data());

            const auto fa = flag4_mf.const_array(mfi);
            auto       oa = out_mf.array(mfi);
            int cnt = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk)
            for (int jj = lo[1]; jj <= hi[1]; ++jj)
            for (int ii = lo[0]; ii <= hi[0]; ++ii)
            {
                if (fa(ii,jj,kk) >= AIR_FLAG) oa(ii,jj,kk) = vals[cnt];
                ++cnt;
            }
        }
    }
#else
    std::vector<double> vals(p->knox * p->knoy * p->knoz);
    count = 0;
    KJILOOP { PFLUIDCHECK vals[count] = in(i,j,k); SFLUIDCHECK vals[count] = 0.0; ++count; }
    HYPRE_SStructVectorSetBoxValues(xv, 0, ilower, iupper, variable, vals.data());
    HYPRE_SStructVectorAssemble(xv);
    HYPRE_SStructVectorAssemble(yv);

    HYPRE_SStructMatrixMatvec(1.0, A, xv, 0.0, yv);

    HYPRE_SStructVectorGetBoxValues(yv, 0, ilower, iupper, variable, vals.data());
    count = 0;
    KJILOOP { PFLUIDCHECK out(i,j,k) = vals[count]; ++count; }
#endif

    HYPRE_SStructVectorDestroy(xv);
    HYPRE_SStructVectorDestroy(yv);
}
