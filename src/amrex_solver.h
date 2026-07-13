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

#if USE_AMREX
#ifndef AMREX_SOLVER_H_
#define AMREX_SOLVER_H_

#include "increment.h"

#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

#include <memory>

class lexer;
class fdm;
class ghostcell;
class field1;
class field2;
class field3;
class field4;
class density;

/// Composite (all-levels-at-once) MAC-projection Poisson solver backed by
/// AMReX MLABecLaplacian/MLMG. setup() stages the face coefficients and the
/// MAC predictor velocity into genuinely face-centred MultiFabs and builds
/// the linear operator; solve() runs the multigrid; ucorr() applies the
/// operator-consistent velocity correction back into the REEF3D fields.
/// start() ties the steps together, mirroring pjm::start.
class amrex_solver : public increment
{
public:
    amrex_solver(lexer *p);
    virtual ~amrex_solver();

    /// Full projection step: setup + rhs + solve + velocity correction +
    /// pressure update. u/v/w and phi need current ghost cells on entry.
    void start(lexer *p, fdm *a, ghostcell *pgc, field1 &u, field2 &v, field3 &w, const field4 &phi, double alpha);

    void setup(lexer *p, fdm *a, ghostcell *pgc, const field1 &u, const field2 &v, const field3 &w, const field4 &phi);

    /// Solves -div(1/ro grad p) = -div(umac)/(alpha*dt) for the staged umac;
    /// returns the final relative residual and sets p->solveriter.
    double solve(lexer *p);

    /// Applies u += alpha*dt*flux with flux = -(1/ro)grad(pcorr) taken from the
    /// operator's own discrete gradient (C-F consistent), refreshes ghost
    /// cells and averages covered coarse faces down.
    void ucorr(lexer *p, fdm *a, ghostcell *pgc, field1 &u, field2 &v, field3 &w, double alpha);

    /// Copies the solved potential into a->press on all levels (the solve is
    /// non-incremental, matching pjm) and refreshes the pressure ghosts.
    void pressure_update(lexer *p, fdm *a, ghostcell *pgc);

private:
    /// rhs = -div(umac)/(alpha*dt) per level from plain face differences.
    void fill_rhs(lexer *p, double alpha);

    density *pd;

    int gcval_press, gcval_u, gcval_v, gcval_w;

    /// Pseudo-2D handling: AMReX's hidden-dimension multi-level path is broken
    /// (divergent composite cycles), so for pseudo-2D runs the solver builds
    /// its own hierarchy in which fine level L has 2^L identical y-planes and
    /// isotropic (2,2,2) refinement -- only standard, well-tested AMReX code
    /// paths are exercised. Level 0 stays ny=1 (its MG chain degenerates to a
    /// pure bottom solve, which is the configuration proven in single-level
    /// runs). Data is staged into plane 0 and replicated; results are read
    /// back from plane 0.
    bool ydouble = false;
    amrex::Vector<amrex::Geometry> sgeom;  ///< solver-side geometry per level
    amrex::Vector<amrex::BoxArray> sgrids; ///< solver-side (y-doubled) grids

    // Face-centred staging fabs in the AMReX low-face convention (face f in
    // direction d is the LOW face of cell f). Reallocated when the hierarchy
    // changes; refilled on every setup() call.
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> umac; ///< MAC predictor velocity
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> beta; ///< 1/ro on faces, 0 on solid faces
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> flux; ///< -(1/ro)grad(pcorr) from getFluxes
    amrex::Vector<amrex::MultiFab> rhs;    ///< cell-centred divergence source
    amrex::Vector<amrex::MultiFab> pcorr;  ///< cell-centred solution, 1 ghost (warm start)

    std::unique_ptr<amrex::MLABecLaplacian> linop;
    std::unique_ptr<amrex::MLMG> mlmg;
};

#endif
#endif
