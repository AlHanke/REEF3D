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

#include"weno_nug_func.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"

void weno_nug_func::iqmin(field& f)
{
    q1 = (f(i-2,j,k)-f(i-3,j,k))/p->DXP[IM3];
    q2 = (f(i-1,j,k)-f(i-2,j,k))/p->DXP[IM2];
    q3 = (f(i,j,k)-f(i-1,j,k))/p->DXP[IM1];
    q4 = (f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
    q5 = (f(i+2,j,k)-f(i+1,j,k))/p->DXP[IP1];
}

void weno_nug_func::jqmin(field& f)
{
    q1 = (f(i,j-2,k)-f(i,j-3,k))/p->DYP[JM3];
    q2 = (f(i,j-1,k)-f(i,j-2,k))/p->DYP[JM2];
    q3 = (f(i,j,k)  -f(i,j-1,k))/p->DYP[JM1];
    q4 = (f(i,j+1,k)-f(i,j,k)  )/p->DYP[JP];
    q5 = (f(i,j+2,k)-f(i,j+1,k))/p->DYP[JP1];
}

void weno_nug_func::kqmin(field& f)
{
    q1 = (f(i,j,k-2)-f(i,j,k-3))/p->DZP[KM3];
    q2 = (f(i,j,k-1)-f(i,j,k-2))/p->DZP[KM2];
    q3 = (f(i,j,k)-f(i,j,k-1))/p->DZP[KM1];
    q4 = (f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
    q5 = (f(i,j,k+2)-f(i,j,k+1))/p->DZP[KP1];
}

void weno_nug_func::iqmax(field& f)
{
    q1 = (f(i-1,j,k)-f(i-2,j,k))/p->DXP[IM2];
    q2 = (f(i,j,k)-f(i-1,j,k))/p->DXP[IM1];
    q3 = (f(i+1,j,k)-f(i,j,k))/p->DXP[IP];
    q4 = (f(i+2,j,k)-f(i+1,j,k))/p->DXP[IP1];
    q5 = (f(i+3,j,k)-f(i+2,j,k))/p->DXP[IP2];
}

void weno_nug_func::jqmax(field& f)
{
    q1 = (f(i,j-1,k)-f(i,j-2,k))/p->DYP[JM2];
    q2 = (f(i,j,k)-f(i,j-1,k))/p->DYP[JM1];
    q3 = (f(i,j+1,k)-f(i,j,k))/p->DYP[JP];
    q4 = (f(i,j+2,k)-f(i,j+1,k))/p->DYP[JP1];
    q5 = (f(i,j+3,k)-f(i,j+2,k))/p->DYP[JP2];
}

void weno_nug_func::kqmax(field& f)
{
    q1 = (f(i,j,k-1)-f(i,j,k-2))/p->DZP[KM2];
    q2 = (f(i,j,k)-f(i,j,k-1))/p->DZP[KM1];
    q3 = (f(i,j,k+1)-f(i,j,k))/p->DZP[KP];
    q4 = (f(i,j,k+2)-f(i,j,k+1))/p->DZP[KP1];
    q5 = (f(i,j,k+3)-f(i,j,k+2))/p->DZP[KP2];
}


void weno_nug_func::isqmin(slice& f)
{
    q1 = (f(i-2,j)-f(i-3,j))/p->DXP[IM3];
    q2 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
    q3 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
    q4 = (f(i+1,j)-f(i,j))/p->DXP[IP];
    q5 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
}

void weno_nug_func::jsqmin(slice& f)
{
    q1 = (f(i,j-2)-f(i,j-3))/p->DYP[JM3];
    q2 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
    q3 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
    q4 = (f(i,j+1)-f(i,j))/p->DYP[JP];
    q5 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
}

void weno_nug_func::isqmax(slice& f)
{
    q1 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
    q2 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
    q3 = (f(i+1,j)-f(i,j))/p->DXP[IP];
    q4 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
    q5 = (f(i+3,j)-f(i+2,j))/p->DXP[IP2];
}

void weno_nug_func::jsqmax(slice& f)
{
    q1 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
    q2 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
    q3 = (f(i,j+1)-f(i,j))/p->DYP[JP];
    q4 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
    q5 = (f(i,j+3)-f(i,j+2))/p->DYP[JP2];
}