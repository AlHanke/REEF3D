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

#ifndef HYPRE_SSAMG_H_
#define HYPRE_SSAMG_H_

#include "solver.h"
#include "increment.h"
#include "vec.h"
#include <_hypre_utilities.h>
#include <HYPRE_sstruct_ls.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <vector>

class hypre_ssamg final : public solver, public increment
{
public:
    hypre_ssamg(lexer*, fdm*, ghostcell*);
    virtual ~hypre_ssamg();

    void start(lexer*, fdm*, ghostcell*, field&, vec&, int) override final;
    void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int) override final {};
    void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final {};
    void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int) override final {};
    void startM(lexer*, ghostcell*, double*, double*, double*, int) override final {};

    void start_solver5(lexer*, fdm*, ghostcell*, field&, vec&, int);

    void solve(lexer*);

    void fill_matrix4(lexer*, fdm*, ghostcell*, field&);
    void fillbackvec4(lexer*, field&, int);

    void make_grid_7p(lexer*, fdm*, ghostcell*);

    void create_solver(lexer*, ghostcell*);
    void delete_solver(lexer*, ghostcell*);

    void amr_graph_entries(lexer*, ghostcell*);

private:
    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencil;
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_SStructSolver   pcg_solver;
    HYPRE_SStructSolver   ssamg;
    HYPRE_SStructVariable vartypes[1];

    int numparts;
    int dimensions;
    int variable;
    int numvar;
    int object_type;

    int ilower[3], iupper[3];
    std::vector<double> values;
    int num_iterations;
    double final_res_norm;
    int stencil_indices[7];
    int nentries;

    int numiter, count, q;

    static constexpr int stencil_size = 7; // 7-point stencil
};

#endif
