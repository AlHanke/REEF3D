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
    if(p->changed)
    {
        destroy_grid();
        make_grid_7p(p, a, pgc);
    }
    #endif

    #if USE_AMREX
    if(!solver_created || created_nlevs != p->nlevs)
    #else
    if(!solver_created)
    #endif
    {
        delete_solver();
        create_solver(p, pgc);
    }

    fill_matrix4(p, a, pgc, f);

    solve(p);

    fillbackvec4(p, f, var);

    std::fill(a->rhsvec.V.begin(), a->rhsvec.V.end(), 0.0);
}
