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

void hypre_ssamg::solve(lexer *p)
{
    p->solveriter = 0;

    // Multi-level: PCG + BoomerAMG on the assembled ParCSR operator (SSAMG cannot set up
    // on multi-part grids; the operator is symmetric, so PCG). Single level keeps SSAMG.
    if (p->nlevs > 1)
    {
        HYPRE_ParCSRPCGSetup(par_solver, par_A, par_b, par_x);
        HYPRE_ParCSRPCGSolve(par_solver, par_A, par_b, par_x);

        HYPRE_PCGGetNumIterations(par_solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(par_solver, &final_res_norm);

        // object_type==HYPRE_PARCSR: refresh the SStruct vector's structured data from
        // the solved ParVector so fillbackvec4's GetBoxValues sees the solution.
        HYPRE_SStructVectorGather(x);
    }
    // N10==40: standalone SSAMG
    else if (p->N10 == 40)
    {
        HYPRE_SStructSSAMGSetup(ssamg, A, b, x);
        HYPRE_SStructSSAMGSolve(ssamg, A, b, x);

        HYPRE_SStructSSAMGGetNumIterations(ssamg, &num_iterations);
        HYPRE_SStructSSAMGGetFinalRelativeResidualNorm(ssamg, &final_res_norm);
    }
    // N10==41: GMRES + SSAMG preconditioner
    else
    {
        HYPRE_SStructGMRESSetup(gmres_solver, A, b, x);
        HYPRE_SStructGMRESSolve(gmres_solver, A, b, x);

        HYPRE_SStructGMRESGetNumIterations(gmres_solver, &num_iterations);
        HYPRE_SStructGMRESGetFinalRelativeResidualNorm(gmres_solver, &final_res_norm);
    }

    p->solveriter = num_iterations;
    p->final_res  = final_res_norm;
}
