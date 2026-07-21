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

#ifndef SOLVER_H_
#define SOLVER_H_

#include <set>

class lexer;
class fdm;
class fdm_fnpf;
class ghostcell;
class field;
class vec;
class matrix_diag;

using namespace std;

class solver
{

public:
    virtual ~solver() = default;

	virtual void start(lexer*,fdm*, ghostcell*, field&, vec&, int)=0;
    virtual void startf(lexer*, ghostcell*, field&, vec&, matrix_diag&, int)=0;
    virtual void startF(lexer*, ghostcell*, double*, vec&, matrix_diag&, int)=0;
    virtual void startV(lexer*, ghostcell*, double*, vec&, matrix_diag&, int)=0;
    virtual void startM(lexer*, ghostcell*, double*, double*, double*, int)=0;

    // Coarse-fine consistent velocity correction (AMR pressure projection). Default
    // is a no-op; only the multi-level capable solver (hypre_ssamg) overrides it.
    // Applied after the per-level interior velocity correction so the face velocities
    // at coarse-fine interfaces are made consistent with the matrix Laplacian.
    virtual void cf_velocity_correction(lexer*, fdm*, ghostcell*,
                                        field&, field&, field&, field&, double) {}

    // Slave the fine-level normal velocity on coarse-fine faces to the coarse value, so the
    // predictor never injects a C-F leak the projection then cannot fully remove. Run on the
    // PREDICTOR velocity (before the rhs divergence). Default no-op; hypre_ssamg overrides.
    virtual void cf_velocity_fill_from_coarse(lexer*, fdm*, ghostcell*,
                                              field&, field&, field&) {}

    // out = A * in, using the assembled pressure matrix (incl. C-F couplings). Default
    // no-op; hypre_ssamg overrides it. Used by the projection-consistency probe to compare
    // L*pcorr against the discrete divergence of the velocity correction.
    virtual void matvec_into(lexer*, fdm*, ghostcell*, field& /*out*/, field& /*in*/) {}

    // Save (save=true) / restore (save=false) the fine LOW C-F normal-face velocities -- the
    // ghost faces (fiv-e) that cf_velocity_correction writes. A start1/2/3 between save and
    // restore (FillPatchTwoLevels + FillCoarseFineNormalGhost) overwrites them with coarse
    // interpolation, wiping the matrix-consistent correction; restore re-asserts it so the
    // fine cell's divergence stays the adjoint of its matrix row. Default no-op; hypre_ssamg
    // overrides (it owns cf_links).
    virtual void cf_lowface_save_restore(lexer*, field&, field&, field&, bool /*save*/) {}

    struct cf_mask
    {
        int level;
        int i,j,k;
        int axis;

        bool operator<(const cf_mask& o) const
        {
            if (level != o.level) return level < o.level;
            if (i     != o.i)     return i     < o.i;
            if (j     != o.j)     return j     < o.j;
            if (k     != o.k)     return k     < o.k;
            return axis < o.axis;
        }
    };

    std::set<cf_mask> cf_masks; // C-F high-faces skipped by the interior velocity loop (velcorr_amrex)
};

#endif
