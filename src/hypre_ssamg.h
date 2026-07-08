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

class fieldint4;

#include <_hypre_utilities.h>
#include <HYPRE_sstruct_ls.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <vector>
#include <map>
#include <array>

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

    void start_solver123(lexer*, fdm*, ghostcell*, field&, int);
    void start_solver45(lexer*, fdm*, ghostcell*, field&, int);

    void solve(lexer*);

    void fill_matrix4(lexer*, fdm*, ghostcell*, field&);
    void amr_cf_coefficients(lexer*, fdm*, ghostcell*, fieldint4&);
    void fillbackvec4(lexer*, field&, int);

    // Momentum diffusion (implicit Helmholtz) solve on the staggered velocity systems.
    // var 1/2/3 -> u/v/w. Single-level only: the velocity fields are cell-centred MultiFabs
    // sharing pressure's 7-point grid, so make_grid_7p is reused. The matrix/RHS come straight
    // from a->M/a->rhsvec as written by idiff2_FS* (single ULOOP/VLOOP/WLOOP pass over
    // flag1/2/3>0), so the numbering is rebuilt with the SAME loop into a sentinel cval field.
    // No C-F coupling / no nullspace projection (Helmholtz is non-singular, diagonally dominant).
    void fill_matrix_vel(lexer*, fdm*, ghostcell*, field&, int var);
    void fillbackvec_vel(lexer*, field&, int var);

    // Sentinel row numbering for the velocity systems: -1 on un-numbered cells, else the
    // sequential index matching idiff2_FS*'s ULOOP(var1)/VLOOP(var2)/WLOOP(var3) fill of
    // a->M/a->rhsvec. Rebuilt identically in fill and fillback so the solved-cell mask agrees.
    void number_velocity(lexer*, fieldint4&, int var);

    void cf_velocity_correction(lexer*, fdm*, ghostcell*,
                                field&, field&, field&, field&, double) override;

    void cf_velocity_fill_from_coarse(lexer*, fdm*, ghostcell*,
                                      field&, field&, field&) override;

    // Operator pre-checks (env REEF_MAT_CHECK): apply the assembled matrix A to test
    // vectors and report on the result. (1) rowsum A*1: ~0 on clean interior rows, the
    // retained BC coefficient on boundary rows. (2) symmetry: global yT A x - xT A y, plus
    // the per-row Ax - A^T x localisation on the multi-level (ParCSR) path.
    void validate_operator(lexer*, fdm*, ghostcell*);

    // out = A * in on the assembled operator (handles SStruct and ParCSR object types).
    void matvec_into(lexer*, fdm*, ghostcell*, field& out, field& in) override;

    // Save/restore the fine LOW C-F normal-face velocities around a start1/2/3 (see solver.h).
    void cf_lowface_save_restore(lexer*, field&, field&, field&, bool save) override;

    void make_grid_7p(lexer*, fdm*, ghostcell*);
    void destroy_grid();

    void create_solver(lexer*, ghostcell*);
    void delete_solver(lexer*, ghostcell*);

private:
    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencil;
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_SStructSolver   gmres_solver;
    HYPRE_SStructSolver   ssamg;
    HYPRE_SStructVariable vartypes[1];

    // Multi-level (nlevs>1) path: the SStruct matrix is assembled as ParCSR and solved
    // with PCG + BoomerAMG (SSAMG cannot set up on multi-part grids). The volume-weighted
    // C-F coupling makes the operator symmetric, so PCG (not BiCGSTAB); the singular
    // all-Neumann system is handled by RHS projection + a non-singular AMG coarse solve.
    HYPRE_ParCSRMatrix    par_A;
    HYPRE_ParVector       par_b;
    HYPRE_ParVector       par_x;
    HYPRE_Solver          par_solver;
    HYPRE_Solver          par_precond;

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

    int count, q;

    static constexpr int stencil_size = 7; // 7-point stencil

    // Coarse-fine couplings (non-stencil graph entries). Built once per regrid in
    // make_grid_7p; .coeff is refreshed every solve during matrix preparation. The
    // .entry index ties each coupling to its HYPRE non-stencil slot (7 + k, where k
    // is the order the entry was added to the graph for the "from" cell), so the
    // matrix fill must replay this list to stay consistent with the graph.
    struct cf_link
    {
        int    from_part;
        int    from_ijk[3];
        int    to_part;
        int    to_ijk[3];
        int    axis;   // 0:x,1:y,2:z
        bool   high;   // low side (s,e,b) or high side (n,w,t) of the from-cell
        int    entry;
        double coeff;
    };
    std::vector<cf_link> cf_links;

    // Scratch for cf_lowface_save_restore: (lev, face_i, face_j, face_k, axis) -> saved velocity.
    std::map<std::array<int,5>,double> cf_lowface_store;
};

#endif
