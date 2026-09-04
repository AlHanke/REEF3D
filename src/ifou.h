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

#ifndef IFOU_H_
#define IFOU_H_

#include"convection.h"
#include"increment.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"
#include<variant>

class ifou final : public convection,  public increment
{
public:
    ifou(lexer*);
    virtual ~ifou() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    template<typename FluxT, typename GenericField>
    void aij(FluxT&, lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&, double*,double*,double*);

    int count;

    std::variant<flux_face_CDS2, flux_face_FOU,
                 flux_face_CDS2_vrans, flux_face_FOU_vrans,
                 flux_face_CDS2_2D, flux_face_FOU_2D,
                 flux_face_CDS2_vrans_2D, flux_face_FOU_vrans_2D> pflux;
};

#endif
