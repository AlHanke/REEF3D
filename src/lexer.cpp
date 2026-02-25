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

#include"lexer.h"

lexer::lexer() : position(this), interpolation(this),
                 flag1(this), flag2(this), flag3(this), flag4(this), flag5(this),
                 flagsf1(this), flagsf2(this), flagsf3(this), flagsf4(this),
                 DF(this), DF1(this), DF2(this), DF3(this),
                 cmu(0.09), sigT(0.9), mpirank(0)
{
    ini_default(this);
}

lexer::~lexer()
{
}
