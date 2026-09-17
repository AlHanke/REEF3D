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

#ifndef POISSON_DUMP_H_
#define POISSON_DUMP_H_

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class field;

// Per-cell dump of the pressure-correction solve, for comparing the AMReX and non-AMReX
// builds cell by cell. Identical output format in both builds; every level and every cell
// is written, INCLUDING covered coarse cells (they are a known defect site, so the mask is
// emitted as a column rather than used as a filter).
//
// Enabled by env REEF_POISSON_DUMP=<tag>. Optional:
//   REEF_POISSON_DUMP_STEP=<n>   only dump step n (default: every step)
//   REEF_POISSON_DUMP_STAGE=<n>  only dump RK stage n, 0-based (default: every stage)
//
// Writes  pdump_<tag>_s<step>_g<stage>_l<lev>_r<rank>.dat  in the run directory, plus
// pdump_<tag>_meta.txt from rank 0. Cell-centre COORDINATES are the join key -- not
// indices, which differ between builds through pseudo-2D y-doubling, per-level offsets,
// margins and the MPI decomposition.
// Derives from increment for the shared static loop indices i,j,k that the LOOP macros
// drive; the coordinate accessors pos_x/y/z live on lexer (via position) and are level-aware.
class poisson_dump : public increment
{
public:
    static void start(lexer *p, fdm *a, ghostcell *pgc, field &pcorr, double alpha);
};

#endif
