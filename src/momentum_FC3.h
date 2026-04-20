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

#ifndef MOMENTUM_FC3_H_
#define MOMENTUM_FC3_H_

#include"momentum.h"
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;
class fluid_update;
class reini;
class picard;
class heat;
class concentration;
class sixdof;
class fsi;

using namespace std;

class momentum_FC3 : public momentum, public momentum_forcing, public bcmom
{
public:
    momentum_FC3(lexer*, fdm*, ghostcell*, convection*, convection*, diffusion*, pressure*, poisson*,
                turbulence*, solver*, solver*, ioflow*, heat*&, concentration*&, reini*, fsi*);
    virtual ~momentum_FC3();
    void start(lexer*, fdm*, ghostcell*, vrans*, sixdof*) override final;

private:
    void irhs(lexer*,fdm*);
    void jrhs(lexer*,fdm*);
    void krhs(lexer*,fdm*);

    #if USE_AMREX
    // Shared MultiFabs for the RK velocity temp fields, declared before the
    // field members so they are constructed first (declaration order).
    // Component layout: 0:u-type  1:v-type  2:w-type
    amrex::Vector<amrex::MultiFab> m_rk1; ///< urk1(0) vrk1(1) wrk1(2)
    amrex::Vector<amrex::MultiFab> m_rk2; ///< urk2(0) vrk2(1) wrk2(2)
    amrex::Vector<amrex::MultiFab> m_f; ///< fx(0) fy(1) fz(2)
    #endif
    field1 udiff,urk1,urk2,fx;
    field2 vdiff,vrk1,vrk2,fy;
    field3 wdiff,wrk1,wrk2,fz;
    field4 ls,frk1,frk2;

    convection *pconvec;
    convection *pfsfdisc;
    diffusion *pdiff;
    pressure *ppress;
    poisson *ppois;
    turbulence *pturb;
    solver *psolv;
    solver *ppoissonsolv;
    ioflow *pflow;
    reini *preini;
    sixdof *p6dof;
    fsi *pfsi;
    fluid_update *pupdate;
    picard *ppicard;

    static constexpr int gcval_u=10, gcval_v=11, gcval_w=12;
    int gcval_phi;
};

#endif
