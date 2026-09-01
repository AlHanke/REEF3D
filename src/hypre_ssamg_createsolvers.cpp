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

    // weighted Jacobi — symmetric, required when used as PCG preconditioner
    HYPRE_SStructSSAMGSetRelaxType(ssamg, 1);
    HYPRE_SStructSSAMGSetRelaxWeight(ssamg, 0.7);

    // V(1,1) cycle
    HYPRE_SStructSSAMGSetNumPreRelax(ssamg, 1);
    HYPRE_SStructSSAMGSetNumPostRelax(ssamg, 1);
    HYPRE_SStructSSAMGSetNumCoarseRelax(ssamg, 2);

    // skip redundant relaxation sweeps on isotropic problems
    HYPRE_SStructSSAMGSetSkipRelax(ssamg, 1);

    // BoomerAMG closes the coarse-level problem
    HYPRE_SStructSSAMGSetCoarseSolverType(ssamg, 1);

    // Hand off to that BoomerAMG coarse solver early instead of driving SSAMG's own
    // structured coarsening to hypre's (much smaller) default. Each extra structured level
    // costs a Galerkin RAP in setup and a relax+restrict+interpolate in every cycle, and on
    // this operator it buys almost no iteration reduction -- BoomerAMG absorbs the same work
    // more cheaply. Measured pressure-solve time on the 2D dam break, iteration count flat:
    //   14.4k cells:  27.1s -> 19.9s   (coarse size 1000)
    //   115.2k cells: 74.0s -> 63.4s   (coarse size 4000)
    // Scaled with the problem so the handoff point stays at the same relative depth; the
    // bounds keep tiny grids from skipping SSAMG entirely (hypre faults with no coarsening
    // at all) and cap the size of the problem BoomerAMG is handed.
    int coarse_size = 500;

    if(p->cellnumtot > 0)
        coarse_size = std::min(std::max(p->cellnumtot / 20, 500), 20000);

    #if USE_AMREX
    // On a multi-part (AMR) grid the handoff has to come much sooner. The C-F couplings live
    // in the graph, not the stencil, so SSAMG's structured coarsening cannot see them and the
    // coarse operators it builds represent the interface badly -- every structured level past
    // the first costs a Galerkin RAP and a full cycle sweep while barely reducing iterations.
    // Handing over after roughly one coarsening (~half the unknowns; cellnumtot counts level 0
    // only) measured 17.3 -> 13.0 GMRES iterations and 52.7s -> 46.0s of pressure-solve time
    // on the 2-level dam break. Note this path is only reachable when the ParCSR branch above
    // is disabled: with it enabled, nlevs>1 never reaches SSAMG, and the ParCSR route is much
    // the faster of the two anyway (8.6 iterations, 29.2s on the same case and build).
    // (nlevs comes from grid_amrex, which only replaces grid in lexer under USE_AMREX; a
    // non-AMReX build is always single level and keeps the rule above.)
    if(p->cellnumtot > 0 && p->nlevs > 1)
        coarse_size = std::min(std::max(p->cellnumtot * 7 / 10, 500), 200000);
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
        HYPRE_SStructGMRESDestroy(gmres_solver);

    HYPRE_SStructSSAMGDestroy(ssamg);
}
