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

#ifndef PJM_CORR_H_
#define PJM_CORR_H_

#include "pressure.h"
#include "pressure_reference.h"
#include "field4.h"

class heat;
class concentration;
class density;
class solver;

using namespace std;

class pjm_corr final : public pressure, public pressure_reference
{
public:
    pjm_corr(lexer*,fdm*,ghostcell*,heat*&,concentration*&);
    virtual ~pjm_corr();

    void start(fdm*,lexer*,poisson*,solver*,ghostcell*,ioflow*,field&,field&,field&,double) override final;
    void ini(lexer*,fdm*,ghostcell*) override final;
    void upgrad(lexer*,fdm*,slice&,slice&) override final;
    void vpgrad(lexer*,fdm*,slice&,slice&) override final;
    void wpgrad(lexer*,fdm*,slice&,slice&) override final;

protected:
    void ucorr(lexer*,fdm*,field&,double) override final;
    void vcorr(lexer*,fdm*,field&,double) override final;
    void wcorr(lexer*,fdm*,field&,double) override final;

private:
    void velcorr(lexer*, fdm*, ghostcell*, field&, field&, field&, solver*, double);
    void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
    void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
    void presscorr(lexer*,fdm*);

    // Projection-consistency probe (env REEF_PROJ_CHECK): applies the full velocity
    // correction to a zero base, then checks L*pcorr + R(dU) ~ 0 per cell. Non-zero
    // residual marks where the discrete divergence of the correction disagrees with
    // the matrix row (the source of the growing spurious velocity).
    void projection_consistency_check(lexer*,fdm*,ghostcell*,solver*,double);
    field4 pcorr;

    density *pd;
};

#endif
