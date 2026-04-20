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

#ifndef MOMENTUM_RKLS3_DF_H_
#define MOMENTUM_RKLS3_DF_H_

#include"momentum.h"
#include"bcmom.h"
#include"diffusion.h"
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
class density;
class poisson;
class sixdof;
class fsi;

using namespace std;

class momentum_RKLS3_df final : public momentum, public bcmom
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    momentum_RKLS3_df(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*, ioflow*);
    virtual ~momentum_RKLS3_df() = default;
    void start(lexer*, fdm*, ghostcell*, vrans*, sixdof*) override final {};

    void starti(lexer*, fdm*, ghostcell*, sixdof*, vrans*, fsi*);

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

    convection *pconvec;
    diffusion *pdiff;
    pressure *ppress;
    poisson *ppois;
    turbulence *pturb;
    solver *psolv;
    solver *ppoissonsolv;
    ioflow *pflow;

    Eigen::Vector3d alpha, gamma, zeta;

    static constexpr int gcval_u=10, gcval_v=11, gcval_w=12;
};

#endif
