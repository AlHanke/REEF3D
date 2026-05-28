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

#ifndef PRESSURE_REFERENCE_F_H_
#define PRESSURE_REFERENCE_F_H_

#include "increment.h"

class lexer;
class fdm;
class ghostcell;

class pressure_reference : virtual public increment
{
public:
    pressure_reference(lexer*);
    virtual ~pressure_reference() = default;

protected:
    void reference_start(lexer*,fdm*,ghostcell*);
    void reference_ini(lexer*,fdm*,ghostcell*);

private:
    void gage_fixed(lexer*,fdm*,ghostcell*);
    void gage_fsf(lexer*,fdm*,ghostcell*);
    void fsf_normalize(lexer*,fdm*,ghostcell*);
};

#endif
