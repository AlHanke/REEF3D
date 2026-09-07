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
    #if USE_AMREX
    if (p->nlevs > 1)
    {
        HYPRE_BoomerAMGCreate(&par_precond);
        HYPRE_BoomerAMGSetPrintLevel(par_precond, 0);
        HYPRE_BoomerAMGSetCoarsenType(par_precond, 22);
        HYPRE_BoomerAMGSetRelaxType(par_precond, 6);     // symmetric hybrid GS
        HYPRE_BoomerAMGSetNumSweeps(par_precond, 1);
        // Coarsen all the way down (9 rows) rather than stopping at 200. Stopping early was
        // only needed to keep BoomerAMG's default Gaussian-elimination coarse solver off a
        // singular grid -- but CycleRelaxType(...,3) below already replaces GE with relaxation,
        // so the early stop bought nothing and left a 200-row coarse problem that one relax
        // sweep cannot solve. Measured on the 2D dam break (2 levels, 18k unknowns): 14.7 -> 8.0
        // GMRES iterations per solve.
        HYPRE_BoomerAMGSetMaxCoarseSize(par_precond, 9);
        HYPRE_BoomerAMGSetCycleRelaxType(par_precond, 6, 3); // relax (not GE) on the coarsest level
        HYPRE_BoomerAMGSetTol(par_precond, 0.0);
        HYPRE_BoomerAMGSetMaxIter(par_precond, 1);

        HYPRE_ParCSRGMRESCreate(pgc->mpi_comm, &par_solver);
        HYPRE_GMRESSetMaxIter(par_solver, p->N46);
        HYPRE_GMRESSetKDim(par_solver, 30);              // restart dimension
        HYPRE_GMRESSetTol(par_solver, p->N44);
        HYPRE_GMRESSetAbsoluteTol(par_solver, 1e-12);
        HYPRE_GMRESSetPrintLevel(par_solver, 0);
        HYPRE_GMRESSetLogging(par_solver, 1);
        HYPRE_GMRESSetPrecond(par_solver,
            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup,
            par_precond);

        // Record the creation state on the multi-level path too. Without this the early
        // return leaves created_nlevs at its stale value (-1) and solver_created false, so
        // solve() takes the created_nlevs<=1 branch and calls HYPRE_SStructGMRESSetup on the
        // never-created single-level gmres_solver (nullptr) -> segfault. delete_solver() also
        // keys off created_nlevs>1 to free the right objects.
        solver_created = true;
        created_nlevs  = p->nlevs;
        grid_rebuilt   = false;

        // This solver has no hierarchy yet, so the next solve must build one before it can
        // start lagging the setup again.
        par_setup_count = 0;
        return;
    }
    #endif

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
    HYPRE_SStructSSAMGSetRelaxType(ssamg, 2);
    HYPRE_SStructSSAMGSetRelaxWeight(ssamg, 1.5);

    // V(1,1) cycle
    HYPRE_SStructSSAMGSetNumPreRelax(ssamg, 1);
    HYPRE_SStructSSAMGSetNumPostRelax(ssamg, 1);
    HYPRE_SStructSSAMGSetNumCoarseRelax(ssamg, 2);

    // skip redundant relaxation sweeps on isotropic problems
    HYPRE_SStructSSAMGSetSkipRelax(ssamg, 1);

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
    HYPRE_SStructSSAMGSetMaxLevels(ssamg, 7);

    // Disable the size-based cutoff so the level cap above is the sole transition criterion,
    // as in the paper. This replaces a locally tuned heuristic (coarse size cellnumtot/20,
    // clamped to [500, 20000], and cellnumtot*7/10 clamped to [500, 200000] for nlevs>1)
    // that handed off far earlier -- ~5x DOF reduction rather than 64x. That heuristic was
    // measured faster on the 2D dam break at the time (14.4k cells 27.1s -> 19.9s, 115.2k
    // cells 74.0s -> 63.4s, iteration count flat), so it is entirely possible SSAMG-opt's
    // deeper structured hierarchy is slower on this operator than what it replaces; the two
    // want to be benchmarked against each other rather than assumed.
    //
    // Note the coarsest grid still ends up at most 9 rows as the paper describes, because
    // ssamg_csolver.c leaves BoomerAMG's own MaxCoarseSize at hypre's default of 9 and does
    // not expose it.
    HYPRE_SStructSSAMGSetMaxCoarseSize(ssamg, 0);

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
    #if USE_AMREX
    created_nlevs = p->nlevs;
    grid_rebuilt  = false;
    #endif
}

void hypre_ssamg::delete_solver()
{
    #if USE_AMREX
    if (created_nlevs > 1)
    {
        HYPRE_ParCSRGMRESDestroy(par_solver);
        HYPRE_BoomerAMGDestroy(par_precond);
        return;
    }
    #endif

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
