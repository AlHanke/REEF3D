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

#include "hypre_ssamg.h"
#include "lexer.h"
#include "ghostcell.h"
#include <algorithm>
void hypre_ssamg::create_solver(lexer *p, ghostcell *pgc)
{
    // ---- Multi-level: ParCSR GMRES + BoomerAMG -----------------------------------
    // SSAMG cannot set up on a multi-part grid, so for nlevs>1 the matrix is assembled
    // as ParCSR (see make_grid_7p) and solved with GMRES preconditioned by one BoomerAMG
    // V-cycle. GMRES (not PCG) because the all-Neumann operator is singular (constant
    // nullspace) and, on a thin adaptive interface band anchored only by the free-surface
    // Dirichlet line, becomes near-singular/ill-conditioned -- PCG+BoomerAMG diverged there
    // (pres ~1e10). GMRES needs no SPD preconditioner and tolerates the near-singular
    // operator via its Krylov least-squares minimisation. Two BoomerAMG cautions remain:
    //   * BoomerAMG's default coarsest-grid solver is Gaussian elimination, which is
    //     singular on the all-Neumann coarse grid and returns garbage. Stop coarsening
    //     early (MaxCoarseSize) and relax the coarsest grid instead of direct-solving it
    //     (CycleRelaxType ..., 3). The RHS is projected onto the compatible subspace in
    //     fill_matrix4, so a min-norm solution exists.
    //   * the symmetric smoother (RelaxType 6) is retained -- harmless for GMRES and keeps
    //     the V-cycle well behaved.
    // par_A/par_b/par_x are extracted after assembly in fill_matrix4; the solver/precond
    // are set up against them in solve().
    // #if USE_AMREX
    // if (p->nlevs > 1)
    // {
    //     HYPRE_BoomerAMGCreate(&par_precond);
    //     HYPRE_BoomerAMGSetPrintLevel(par_precond, 0);
    //     HYPRE_BoomerAMGSetCoarsenType(par_precond, 22);
    //     HYPRE_BoomerAMGSetRelaxType(par_precond, 6);     // symmetric hybrid GS
    //     HYPRE_BoomerAMGSetNumSweeps(par_precond, 1);
    //     // Coarsen all the way down (9 rows) rather than stopping at 200. Stopping early was
    //     // only needed to keep BoomerAMG's default Gaussian-elimination coarse solver off a
    //     // singular grid -- but CycleRelaxType(...,3) below already replaces GE with relaxation,
    //     // so the early stop bought nothing and left a 200-row coarse problem that one relax
    //     // sweep cannot solve. Measured on the 2D dam break (2 levels, 18k unknowns): 14.7 -> 8.0
    //     // GMRES iterations per solve.
    //     HYPRE_BoomerAMGSetMaxCoarseSize(par_precond, 9);
    //     HYPRE_BoomerAMGSetCycleRelaxType(par_precond, 6, 3); // relax (not GE) on the coarsest level
    //     HYPRE_BoomerAMGSetTol(par_precond, 0.0);
    //     HYPRE_BoomerAMGSetMaxIter(par_precond, 1);

    //     HYPRE_ParCSRGMRESCreate(pgc->mpi_comm, &par_solver);
    //     HYPRE_GMRESSetMaxIter(par_solver, p->N46);
    //     HYPRE_GMRESSetKDim(par_solver, 30);              // restart dimension
    //     HYPRE_GMRESSetTol(par_solver, p->N44);
    //     HYPRE_GMRESSetAbsoluteTol(par_solver, 1e-12);
    //     HYPRE_GMRESSetPrintLevel(par_solver, 0);
    //     HYPRE_GMRESSetLogging(par_solver, 1);
    //     HYPRE_GMRESSetPrecond(par_solver,
    //         (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
    //         (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
    //         par_precond);

    //     // Record the creation state on the multi-level path too. Without this the early
    //     // return leaves created_nlevs at its stale value (-1) and solver_created false, so
    //     // solve() takes the created_nlevs<=1 branch and calls HYPRE_SStructGMRESSetup on the
    //     // never-created single-level gmres_solver (nullptr) -> segfault. delete_solver() also
    //     // keys off created_nlevs>1 to free the right objects.
    //     solver_created = true;
    //     created_nlevs  = p->nlevs;
    //     grid_rebuilt   = false;

    //     // This solver has no hierarchy yet, so the next solve must build one before it can
    //     // start lagging the setup again.
    //     setup_count = 0;
    //     return;
    // }
    // #endif

    // ---- Single level: SSAMG (native SStruct) -----------------------------------
    // SSAMG preconditioner / standalone solver
    HYPRE_SStructSSAMGCreate(pgc->mpi_comm, &ssamg);

    // structured-only interpolation within parts (fastest for block-structured grids)
    HYPRE_SStructSSAMGSetInterpType(ssamg, -1);

    // Smoother: weighted L1-Jacobi with relaxation factor 3/2 -- the SSAMG-opt setting of
    // Magri, Falgout & Yang, "A New Semistructured Algebraic Multigrid Method",
    // SIAM J. Sci. Comput. 45(3), S439-S460 (2023), sec. 5. Relax type 2 is L1-Jacobi
    // (ssamg_relax.c builds z from hypre_SStructMatrixComputeL1Norms instead of the plain
    // diagonal). L1-Jacobi guarantees rho(I - M^-1 A) < 1 for SPD A, and the paper notes a
    // user weight in (1, 2/lambda_max(M^-1 A)) -- of which 3/2 is their choice -- recovers
    // the convergence L1-Jacobi otherwise gives up against weighted Jacobi. Still a diagonal
    // smoother, so it stays symmetric and remains valid as a PCG preconditioner.
    //
    // Caveat worth knowing: calling SetRelaxWeight at all sets usr_set_rweight, which makes
    // ssamg_setup.c:869 skip hypre's automatic per-level/per-part weight
    // omega_p = 2/(3 - beta_p/alpha_p) (the paper's eq. 4.2, which adapts the weight to the
    // anisotropy of each part) and pin this single value on every level and part instead.
    // That is what the paper does for its L1-Jacobi variants; the adaptive formula targets
    // plain weighted Jacobi, where it yields values below one.
    //
    // Swept on the wave-over-bar case and this pair won, so the paper's choice carries over to
    // this operator. L1-Jacobi at w = 1.0 / 1.25 / 1.5 / 1.75 gave 20.7 / 19.8 / 19.1 / 19.3 s,
    // and letting hypre pick the weight automatically (by not calling SetRelaxWeight) gave
    // 20.6 s. Plain weighted Jacobi (relax type 1) was never better -- 19.5 s at w = 0.7,
    // 19.1 s at w = 0.85 -- and at w = 1.0 it DIVERGES: 60 iterations per solve against a
    // 250 cap, a final residual of 2.7, and a visibly wrong free surface. Weighted Jacobi has
    // no damping margin left at w = 1 on this operator, which is exactly the failure L1-Jacobi
    // is chosen to avoid.
    HYPRE_SStructSSAMGSetRelaxType(ssamg, 2);
    HYPRE_SStructSSAMGSetRelaxWeight(ssamg, 1.5);

    // V(1,1) cycle. Swept: (1,1) 19.1s, (2,2) 19.4s, (2,1) 20.5s, (1,2) 20.9s. The heavier
    // cycles do cut iterations -- (2,2) reaches 8.3 against (1,1)'s 10.3 -- but not by enough
    // to pay for the extra sweeps, and the asymmetric pairs lose on both counts.
    HYPRE_SStructSSAMGSetNumPreRelax(ssamg, 1);
    HYPRE_SStructSSAMGSetNumPostRelax(ssamg, 1);

    // One sweep on the coarsest structured level, not two. With BoomerAMG closing the coarse
    // problem below it (CoarseSolverType 1), a second sweep here is redundant work on the
    // level BoomerAMG is about to solve properly: 1 sweep 18.3s, 2 sweeps 19.1s, 4 sweeps
    // 21.6s, 8 sweeps 26.2s, while the iteration count barely moves (11.3 / 10.2 / 10.2).
    HYPRE_SStructSSAMGSetNumCoarseRelax(ssamg, 1);

    // SkipRelax is hypre's isotropy shortcut: ssamg_setup.c marks a level INACTIVE unless its
    // coarsening direction has already been coarsened once before, so relaxation happens about
    // once per full sweep of all directions instead of after every semicoarsening step. How
    // much that actually skips depends on how many directions coarsen -- one level in two on a
    // pseudo-2D grid, two in three on a full 3D one -- so it is not worth the same thing in the
    // two cases, and the sign of the effect genuinely flips. Hence the key on j_dir, the same
    // pseudo-2D discriminator make_grid_7p uses.
    //
    // Measured both ways on both geometries, min of 3, everything else at the values set here:
    //   pseudo-2D wave-over-bar (1200x1x40, 2 lev):  skip 0 = 18.1s / 7.7 it
    //                                                skip 1 = 19.1s / 10.3 it   -> skip costs 5%
    //   3D wave-over-bar        (600x10x20, 2 lev):  skip 0 = 75.6s / 10.1 it
    //                                                skip 1 = 74.9s / 12.6 it   -> skip gains 1.3%
    // Skipping always raises the iteration count; the question is only whether the sweeps it
    // saves are worth more. In pseudo-2D just x and z coarsen, so skipping halves the smoothing
    // on a grid whose free-surface density jump already makes the smoother work hard, and it
    // does not pay. In 3D the third direction leaves enough smoothing per cycle that it does.
    // Both deltas are small but reproduced across three independent measurements each.
    HYPRE_SStructSSAMGSetSkipRelax(ssamg, p->j_dir ? 1 : 0);

    // BoomerAMG closes the coarse-level problem
    HYPRE_SStructSSAMGSetCoarseSolverType(ssamg, 1);

    // "Hybrid" handoff to that BoomerAMG coarse solver, at the depth the paper calls
    // SSAMG-opt: transition at the 7th level, i.e. six pure SSAMG coarsening levels and a
    // 64x reduction in DOFs before BoomerAMG takes over. SSAMG semicoarsens one direction
    // per level, so six levels is two full coarsenings of a 3D grid (three of a 2D one).
    // SSAMG-opt differs from the paper's SSAMG-hybrid in this number alone -- hybrid hands
    // off at the 10th level (512x) -- and it was the fastest of their four variants.
    //
    // hypre exposes no transition-level setter, so the level cap is the mechanism:
    // ssamg_setup.c stops coarsening at l == max_levels-1 and ssamg_csolver.c then converts
    // that coarsest SStructMatrix to ParCSR and hands it to BoomerAMG. Level 7 is only an
    // upper bound -- the loop still exits early once no part has a coarsenable direction
    // left, so small or thin (pseudo-2D) grids simply get fewer levels.
    //
    // With the size cutoff below restored, that cutoff nearly always bites first and this cap
    // is inert: sweeping it over 3/4/5/7/10 at coarse size 20000 moved the total by less than
    // the run-to-run spread (19.1 / 19.07 / 18.94 / 18.95 / 19.08 s) and did not change the
    // iteration count at all above 3. Kept at the paper's value as a backstop.
    HYPRE_SStructSSAMGSetMaxLevels(ssamg, 7);

    // Hand off to the BoomerAMG coarse solver before SSAMG's own structured coarsening reaches
    // hypre's default. Each extra structured level costs a Galerkin RAP in setup and a
    // relax+restrict+interpolate in every cycle, and past a point that buys fewer iterations
    // than it costs. This restores (and re-tunes) the size-based cutoff that the SSAMG-opt
    // configuration of Magri, Falgout & Yang had disabled with MaxCoarseSize(0); the level cap
    // below is then only a backstop.
    //
    // Measured on two geometries, min of 3, with the setup lag in place. Total wall time and
    // mean PCG iterations against coarse size:
    //   pseudo-2D (48k level-0 cells + 77k on level 1, 100 steps, 4 ranks)
    //        0 (off) 22.8s/15.1   2000 18.8s/10.2   5000 18.2s/9.5
    //       10000    18.5s/10.2  20000 17.8s/8.9   40000 18.5s/9.3
    //   3D (120k level-0 cells, 30 steps, 6 ranks)
    //        0 (off) 73.4s/13.7   1000 73.3s/13.7   3000 73.4s/13.7   6000 73.3s/13.7
    //       24000    73.4s/13.1  48000 74.9s/12.6
    // The two disagree about what the cutoff is worth -- it is a 22% win in pseudo-2D and a
    // wash in 3D, where the curve is flat until it starts costing around 48k -- but they agree
    // on where to put it: ~20k is at the pseudo-2D optimum and still inside the 3D flat region.
    //
    // Hence a fraction of the grid with a hard absolute cap, and the CAP is the part that
    // carries most of the generality. BoomerAMG's cost scales with the absolute size of the
    // problem handed to it, not with what fraction of the original grid that is, so "coarsen
    // until at most ~20k rows remain" is the meaningful rule; the fraction only stops small
    // grids from handing over a problem that is most of the grid. An earlier version of this
    // used the fraction alone (cellnumtot*2/5, no cap), which is fine at 48k but picks 48000 on
    // the 3D case -- precisely the value measured to be the worst of the six tested there.
    //
    // The lower clamp keeps tiny grids from skipping SSAMG altogether (hypre faults when
    // nothing coarsens at all). The upper clamp is not cosmetic either: a MaxCoarseSize at or
    // above the whole grid segfaults in the coarse solve (observed at 200000 on the 48k case,
    // SIGSEGV on rank 2), so it has to stay well clear of the grid size.
    int coarse_size = 500;

    if(p->cellnumtot > 0)
        coarse_size = std::min(std::max(p->cellnumtot * 2 / 5, 500), 20000);

    #if USE_AMREX
    // Single-part grids were tuned separately and earlier, on the 2D dam break, where a much
    // earlier handoff (cellnumtot/20) measured best: 14.4k cells 27.1s -> 19.9s, 115.2k cells
    // 74.0s -> 63.4s. Nothing in this round re-measured the single-level case -- every run was
    // 2-level -- so leave that result standing rather than overwriting it with a number tuned
    // on a different operator.
    if(p->cellnumtot > 0 && p->nlevs == 1)
        coarse_size = std::min(std::max(p->cellnumtot / 20, 500), 20000);
    #endif

    HYPRE_SStructSSAMGSetMaxCoarseSize(ssamg, coarse_size);

    // Galerkin RAP works for both single-level and multi-level grids.
    // Non-Galerkin RAP keeps a compact stencil on coarse levels but crashes
    // with a single SStruct part (no inter-part graph entries): the internal
    // IJMatrix is left uninitialised and assembly faults.  Re-enable once
    // multi-level AMR is confirmed working.
    // if (numparts > 1)
    //     HYPRE_SStructSSAMGSetNonGalerkinRAP(ssamg, 1);
    // else
        HYPRE_SStructSSAMGSetNonGalerkinRAP(ssamg, 0);

    HYPRE_SStructSSAMGSetLogging(ssamg, 0);
    HYPRE_SStructSSAMGSetPrintLevel(ssamg, 0);

    // N10==40: standalone SSAMG direct solver
    if (p->N10 == 40)
    {
        HYPRE_SStructSSAMGSetMaxIter(ssamg, p->N46);
        HYPRE_SStructSSAMGSetTol(ssamg, p->N44);
        HYPRE_SStructSSAMGSetNonZeroGuess(ssamg);
    }
    // N10==41: GMRES outer solver with SSAMG preconditioner (one V-cycle per iteration)
    else if (p->N10 == 41)
    {
        HYPRE_SStructSSAMGSetMaxIter(ssamg, 1);
        HYPRE_SStructSSAMGSetTol(ssamg, 0.0);
        HYPRE_SStructSSAMGSetZeroGuess(ssamg);

        HYPRE_SStructGMRESCreate(pgc->mpi_comm, &gmres_solver);
        HYPRE_SStructGMRESSetMaxIter(gmres_solver, p->N46);
        HYPRE_SStructGMRESSetKDim(gmres_solver, 30);
        HYPRE_SStructGMRESSetTol(gmres_solver, p->N44);
        HYPRE_SStructGMRESSetPrintLevel(gmres_solver, 0);
        HYPRE_SStructGMRESSetLogging(gmres_solver, 1);

        HYPRE_SStructGMRESSetPrecond(gmres_solver,
            HYPRE_SStructSSAMGSolve,
            HYPRE_SStructSSAMGSetup,
            ssamg);

        gmres_created = true;
    }
    // N10==42: PCG outer solver with SSAMG preconditioner -- the setup every result in the
    // SSAMG paper is measured with. SSAMG is never run standalone there; §5 states the
    // preconditioner "is applied to the residual vector via a single V(1,1)-cycle", which is
    // exactly the MaxIter(1)/Tol(0)/ZeroGuess configuration below, and the paper's stopping
    // criterion is ||r||_2 < 1e-6 ||b||_2 from a zero initial guess.
    //
    // WARNING -- read before selecting this. PCG requires BOTH the operator and the
    // preconditioner to be SPD, and neither is guaranteed here:
    //   * The paper's test problem is a Poisson system with Dirichlet data on the k=0 face,
    //     so it is nonsingular and SPD. REEF3D's pressure Poisson is all-Neumann and singular
    //     (constant nullspace) whenever no free surface pins it, and near-singular on a thin
    //     interface band. PCG has no defence against that: the multi-level path above uses
    //     GMRES precisely because PCG+BoomerAMG diverged there (press ~1e10).
    //   * SSAMG stays symmetric only because the smoother is diagonal (L1-Jacobi) and the
    //     cycle is V(1,1) -- so do not switch the smoother to red/black Gauss-Seidel
    //     (relax type 10) on this path without symmetrising the cycle.
    // N10==41 (GMRES) is the safe default for this solver; 42 exists to reproduce the paper's
    // configuration on problems that are actually SPD, and to A/B the SSAMG-opt parameters
    // under the Krylov method they were tuned for.
    else if (p->N10 == 42)
    {
        HYPRE_SStructSSAMGSetMaxIter(ssamg, 1);
        HYPRE_SStructSSAMGSetTol(ssamg, 0.0);
        HYPRE_SStructSSAMGSetZeroGuess(ssamg);

        HYPRE_SStructPCGCreate(pgc->mpi_comm, &pcg_solver);
        HYPRE_SStructPCGSetMaxIter(pcg_solver, p->N46);
        HYPRE_SStructPCGSetTol(pcg_solver, p->N44);
        // Stop on the true residual 2-norm, matching the paper's ||r||_2 < tol*||b||_2
        // rather than PCG's default preconditioned norm.
        HYPRE_SStructPCGSetTwoNorm(pcg_solver, 1);
        HYPRE_SStructPCGSetRelChange(pcg_solver, 0);
        HYPRE_SStructPCGSetPrintLevel(pcg_solver, 0);
        HYPRE_SStructPCGSetLogging(pcg_solver, 1);

        HYPRE_SStructPCGSetPrecond(pcg_solver,
            HYPRE_SStructSSAMGSolve,
            HYPRE_SStructSSAMGSetup,
            ssamg);

        pcg_created = true;
    }

    solver_created = true;
    // Cleared unconditionally: a non-AMReX build still sets grid_rebuilt in make_grid_7p, and
    // leaving it set would force a rebuild on every solve and defeat the lag.
    grid_rebuilt = false;
    #if USE_AMREX
    created_nlevs = p->nlevs;
    #endif
}

void hypre_ssamg::delete_solver()
{
    // #if USE_AMREX
    // if (created_nlevs > 1)
    // {
    //     HYPRE_ParCSRGMRESDestroy(par_solver);
    //     HYPRE_BoomerAMGDestroy(par_precond);
    //     return;
    // }
    // #endif

    if (gmres_created)
    {
        HYPRE_SStructGMRESDestroy(gmres_solver);
        gmres_created = false;
    }

    if (pcg_created)
    {
        HYPRE_SStructPCGDestroy(pcg_solver);
        pcg_created = false;
    }

    HYPRE_SStructSSAMGDestroy(ssamg);
}
