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
#include "fdm.h"
#include "ghostcell.h"
#include "field.h"
#include "vec.h"

hypre_ssamg::hypre_ssamg(lexer *p, fdm *a, ghostcell *pgc)
{
    make_grid_7p(p, a, pgc);
}

hypre_ssamg::~hypre_ssamg()
{
    delete_solver();
    destroy_grid();
}

void hypre_ssamg::start(lexer *p, fdm *a, ghostcell *pgc, field &f, vec&, int var)
{
    if(var==1)
        start_solver123(p, a, pgc, f, var);
    else if(var==2)
        start_solver123(p, a, pgc, f, var);
    else if(var==3)
        start_solver123(p, a, pgc, f, var);
    else if (var == 4 || var == 5)
        start_solver45(p, a, pgc, f, var);
}

void hypre_ssamg::start_solver45(lexer *p, fdm *a, ghostcell *pgc, field &f, int var)
{
    p->solveriter = 0;

    #if USE_AMREX
    // Rebuild on a regrid, and also whenever the assembled grid's part count no longer
    // matches p->nlevs. p->changed is a one-shot flag that lexer::regrid CLEARS on entry:
    // the init hierarchy build (driver_ini_cfd) runs regrid in a loop whose last, no-growth
    // pass resets it to false, so the level 1->2 creation that happened in the previous pass
    // is never seen here -- the solver objects are constructed by logic_cfd() BEFORE that
    // loop, at nlevs==1. fill_matrix4/fill_matrix_vel then loop lev<p->nlevs and call
    // HYPRE_SStructMatrixSetBoxValues on a part that was never given any extents -> segfault.
    if(p->changed || numparts != p->nlevs)
    {
        destroy_grid();
        make_grid_7p(p, a, pgc);
    }
    #endif

    // Rebuild the solver for every solve. The matrix values change every step, so solve() has
    // to run HYPRE_SStructGMRESSetup/SSAMGSetup each time -- and hypre's SStruct setup allocates
    // a complete new hierarchy on every call (SSAMGSetup -> ComputeRAP / MatvecSetup ->
    // hypre_StructMatrixResize) without releasing the one from the previous call. There is no
    // "unsetup", only Destroy, so a solver object kept across steps leaks one full multigrid
    // hierarchy per solve (measured ~25 MB/step at 240x2x240 with 3 RK3 pressure solves, >5 GB
    // after ~200 steps). Create/Destroy are trivial next to the Setup that has to run anyway.
    bool rebuild_solver = true;

    #if USE_AMREX
    // nlevs>1 takes the ParCSR GMRES + BoomerAMG path, which does NOT leak: hypre_BoomerAMGSetup
    // frees its own previous hierarchy. There the rebuild only has to follow the operator's
    // identity -- grid_rebuilt (make_grid_7p destroyed A/b/x, and a solver kept across that still
    // points at freed matrix data -> segfault in hypre_BoomerAMGCycle) or a change in level count.
    if(p->nlevs > 1)
        rebuild_solver = (!solver_created || created_nlevs != p->nlevs || grid_rebuilt);
    #endif

    if(rebuild_solver)
    {
        delete_solver();
        create_solver(p, pgc);
    }

    fill_matrix4(p, a, pgc, f);

    solve(p);

    fillbackvec4(p, f, var);

    std::fill(a->rhsvec.V.begin(), a->rhsvec.V.end(), 0.0);
}
