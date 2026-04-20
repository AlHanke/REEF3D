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
Author: Tobias Martin
--------------------------------------------------------------------*/

#ifndef MOMENTUM_FCLS3_H_
#define MOMENTUM_FCLS3_H_

#include"momentum.h"
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include<vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class fluid_update;
class picard;
class density;
class heat;
class concentration;
class reini;
class poisson;
class sixdof_base;
class fsi;

using namespace std;

class momentum_FCLS3 final : public momentum, public momentum_forcing, public bcmom
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    momentum_FCLS3(lexer*, fdm*, ghostcell*, convection*, convection*, diffusion*, pressure*, poisson*,
                turbulence*, solver*, solver*, ioflow*, heat*&, concentration*&, reini*, fsi*);
    virtual ~momentum_FCLS3();
    void start(lexer*, fdm*, ghostcell*, vrans*, sixdof*) override final;

private:

    void irhs(lexer*,fdm*);
    void jrhs(lexer*,fdm*);
    void krhs(lexer*,fdm*);

    #if USE_AMREX
    // Shared MultiFabs for the RK velocity temp fields, declared before the
    // field members so they are constructed first (declaration order).
    // Component layout: 0:u-type  1:v-type  2:w-type
    amrex::Vector<amrex::MultiFab> m_rk; ///< urk(0) vrk(1) wrk(2)
    amrex::Vector<amrex::MultiFab> m_f; ///< fx(0) fy(1) fz(2)
    #endif
    field1 urk, Cu, fx;
    field2 vrk, Cv, fy;
    field3 wrk, Cw, fz;

    field4 Cf;

    convection *pconvec;
    convection *pfsfdisc;
    diffusion *pdiff;
    diffusion *pdiff_e;
    pressure *ppress;
    poisson *ppois;
    density *pdensity;
    turbulence *pturb;
    solver *psolv;
    solver *ppoissonsolv;
     reini *preini;
    ioflow *pflow;
    fsi *pfsi;
    fluid_update *pupdate;
    picard *ppicard;

    Eigen::Vector3d alpha, gamma, zeta;

    static constexpr int gcval_u=10, gcval_v=11, gcval_w=12;
    int gcval_phi;
};

#endif

