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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef WENO_HJ_DF_NUG_H_
#define WENO_HJ_DF_NUG_H_

#include"convection.h"
#include"weno_nug_func.h"

class flux;

class weno_hj_df_nug final : public convection, public weno_nug_func
{
public:
    weno_hj_df_nug(lexer*);
    virtual ~weno_hj_df_nug() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    template<typename GenericField>
    inline double aij(lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&, double*, double*, double*);

    template<typename GenericField>
    double fx(lexer*, fdm*, const GenericField&, const double*, double);
    template<typename GenericField>
    double fy(lexer*, fdm*, const GenericField&, const double*, double);
    template<typename GenericField>
    double fz(lexer*, fdm*, const GenericField&, const double*, double);

    flux *pflux;
};

#endif
