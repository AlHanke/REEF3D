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
    HYPRE_Int iters; HYPRE_Real relres;

    // Multi-level: GMRES + BoomerAMG on the assembled ParCSR operator (SSAMG cannot set up
    // on multi-part grids; the near-singular all-Neumann operator needs GMRES, not PCG).
    // Single level keeps SSAMG.
    if (created_nlevs > 1)
    {
        // Only rebuild the BoomerAMG hierarchy periodically (see par_setup_count in the
        // header): on a fresh solver, on a fixed period, or as soon as the iteration count
        // shows the lagged hierarchy going stale. GMRES always applies the current par_A, so
        // reusing an older hierarchy changes only how fast it converges, not what it converges
        // to -- the dam break solutions are bit-identical to rebuilding every solve. Setup fell
        // from 13.5s to ~1.0s over 1476 solves for +0.6 iterations each.
        const bool do_setup = (par_setup_count == 0)
                           || (par_setup_count >= par_setup_period)
                           || (par_last_iters > par_fresh_iters + par_setup_degrade);

        if(do_setup)
        {
            HYPRE_ParCSRGMRESSetup(par_solver, par_A, par_b, par_x);
            par_setup_count = 0;
        }
        ++par_setup_count;

        HYPRE_ParCSRGMRESSolve(par_solver, par_A, par_b, par_x);

        HYPRE_GMRESGetNumIterations(par_solver, &iters);
        HYPRE_GMRESGetFinalRelativeResidualNorm(par_solver, &relres);

        par_last_iters = int(iters);
        if(do_setup)
            par_fresh_iters = int(iters);

        // object_type==HYPRE_PARCSR: refresh the SStruct vector's structured data from
        // the solved ParVector so fillbackvec4's GetBoxValues sees the solution.
        HYPRE_SStructVectorGather(x);
    }
    // N10==40: standalone SSAMG
    else if (p->N10 == 40)
    {
        HYPRE_SStructSSAMGSetup(ssamg, A, b, x);
        HYPRE_SStructSSAMGSolve(ssamg, A, b, x);

        HYPRE_SStructSSAMGGetNumIterations(ssamg, &iters);
        HYPRE_SStructSSAMGGetFinalRelativeResidualNorm(ssamg, &relres);
    }
    // N10==41: GMRES + SSAMG preconditioner
    else
    {
        HYPRE_SStructGMRESSetup(gmres_solver, A, b, x);
        HYPRE_SStructGMRESSolve(gmres_solver, A, b, x);

        HYPRE_SStructGMRESGetNumIterations(gmres_solver, &iters);
        HYPRE_SStructGMRESGetFinalRelativeResidualNorm(gmres_solver, &relres);
    }

    p->solveriter = int(iters);
    p->final_res  = double(relres);
}
