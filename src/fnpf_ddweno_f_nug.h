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

#ifndef FNPF_DDWENO_F_NUG_H_
#define FNPF_DDWENO_F_NUG_H_

#include"weno_nug_func.h"

class lexer;
class fdm_fnpf;
class field;
class slice;

using namespace std;

class fnpf_ddweno_f_nug : public weno_nug_func
{
public:
    fnpf_ddweno_f_nug(lexer*,fdm_fnpf*);
    ~fnpf_ddweno_f_nug();

    // field
    double ddwenox(field&, double);
    double ddwenoy(field&, double);
    double ddwenoz(field&, double);

    // slice
    double dswenox(slice&, double);
    double dswenoy(slice&, double);

private:
    inline void iqmin(field&);
    inline void jqmin(field&);
    inline void kqmin(field&);
    inline void iqmax(field&);
    inline void jqmax(field&);
    inline void kqmax(field&);

    inline void isqmin(slice&);
    inline void jsqmin(slice&);
    inline void isqmax(slice&);
    inline void jsqmax(slice&);

    double *DX,*DY,*DZ;

    lexer *p;
    fdm_fnpf *c;
};

#endif
