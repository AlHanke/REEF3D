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

#include"resize.h"

void resize_class::Darray(double*& field, int numi)
{
    if(numi > 0)
    {
        field = new double[numi] {0.0};
    }
}

void resize_class::Darray(double**& field, int numi, int numj)
{    
    if(numi > 0)
    {
        field = new double*[numi];

        if(numj > 0)
            for(int n = 0; n < numi; ++n)
            {
                field[n] = new double[numj] {0.0};
            }
    }
}

void resize_class::Iarray(int*& field, int numi)
{
    if(numi>0)
    {
        field = new int[numi] {0};
    }
}

void resize_class::Iarray(int**& field, int numi, int numj)
{
    if(numi>0)
    {
        field = new int*[numi];

        if(numj > 0)
            for(int n = 0; n < numi; ++n)
            {
                field[n] = new int[numj] {0};
            }
    }
}
