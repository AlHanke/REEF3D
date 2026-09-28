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

#ifndef DENSITY_H_
#define DENSITY_H_

#include"increment.h"

class lexer;
class fdm;

// virtual public increment so update_faces() can use the loop macros; every concrete
// density already inherits increment virtually, so the base stays unambiguous.
class density : virtual public increment
{
public:
    density() = default;
    density(const density&) = delete;
    density& operator=(const density&) = delete;
    density(density&&) = delete;
    density& operator=(density&&) = delete;
    virtual ~density() = default;

    virtual double roface(lexer*,fdm*,int,int,int)=0;

    /// Materialise roface() into a->rofx/rofy/rofz and, with AMR levels, make the
    /// value single-valued across every coarse-fine face. Non-virtual: it drives the
    /// subclass roface(), so all eight density models inherit it unchanged.
    void update_faces(lexer*, fdm*);

    /// Face density from the level-set values at the two ends of an arbitrary segment. Used for
    /// the fine side of a coarse-fine face, whose segment runs from the fine cell centre to the
    /// coarse cell centre (d_cf), not to the interpolated fine ghost. Returns false for models
    /// without such a form; update_faces then keeps roface() there.
    virtual bool roface_segment(lexer*, double, double, double&) { return false; }

private:
    void cf_face_density(lexer*, fdm*);
};

#endif
