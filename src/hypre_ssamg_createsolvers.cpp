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

#include"hypre_ssamg.h"
#include"lexer.h"
#include"ghostcell.h"

void hypre_ssamg::create_solver(lexer *p, ghostcell *pgc)
{
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

    // non-Galerkin RAP keeps compact stencil on coarse levels
    HYPRE_SStructSSAMGSetNonGalerkinRAP(ssamg, 1);

    HYPRE_SStructSSAMGSetLogging(ssamg, 0);
    HYPRE_SStructSSAMGSetPrintLevel(ssamg, 0);

    // N10==40: standalone SSAMG direct solver
    if (p->N10 == 40)
    {
        HYPRE_SStructSSAMGSetMaxIter(ssamg, p->N46);
        HYPRE_SStructSSAMGSetTol(ssamg, p->N44);
        HYPRE_SStructSSAMGSetNonZeroGuess(ssamg);
    }

    // N10==41: PCG outer solver with SSAMG preconditioner (one V-cycle per iteration)
    if (p->N10 == 41)
    {
        HYPRE_SStructSSAMGSetMaxIter(ssamg, 1);
        HYPRE_SStructSSAMGSetTol(ssamg, 0.0);
        HYPRE_SStructSSAMGSetZeroGuess(ssamg);

        HYPRE_SStructPCGCreate(pgc->mpi_comm, &pcg_solver);
        HYPRE_SStructPCGSetMaxIter(pcg_solver, p->N46);
        HYPRE_SStructPCGSetTol(pcg_solver, p->N44);
        HYPRE_SStructPCGSetTwoNorm(pcg_solver, 1);
        HYPRE_SStructPCGSetRelChange(pcg_solver, 0);
        HYPRE_SStructPCGSetLogging(pcg_solver, 1);
        HYPRE_SStructPCGSetPrintLevel(pcg_solver, 0);

        HYPRE_SStructPCGSetPrecond(pcg_solver,
            HYPRE_SStructSSAMGSolve,
            HYPRE_SStructSSAMGSetup,
            ssamg);
    }
}

void hypre_ssamg::delete_solver(lexer *p, ghostcell *pgc)
{
    if (p->N10 == 41)
        HYPRE_SStructPCGDestroy(pcg_solver);

    HYPRE_SStructSSAMGDestroy(ssamg);
}
