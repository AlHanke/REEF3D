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

#ifndef EDIFF2_H_
#define EDIFF2_H_

#include "diffusion.h"
#include "gradient.h"

class ediff2 final : public diffusion, public gradient
{
public:
    ediff2(lexer*);
    virtual ~ediff2() = default;

    void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double, double) override final;
    void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double, double) override final;
    void idiff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, double, double) override final;

    void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double) override final;
    void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double) override final;
    void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double) override final;

private:
    static constexpr int gcval_u = 10, gcval_v = 11, gcval_w = 12;
};
#endif
