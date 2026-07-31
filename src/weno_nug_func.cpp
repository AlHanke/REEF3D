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

weno_nug_func::weno_nug_func(lexer* p):epsilon(0.0),psi(1.0e-6)
{
    ini(p);

    weno_nug_func::p=p;
}

weno_nug_func::~weno_nug_func()
{
}

void weno_nug_func::ini(lexer* p)
{
    if(iniflag==0)
    {
    qfx.resize(p->knox+8);
    qfy.resize(p->knoy+8);
    qfz.resize(p->knoz+8);
    
    cfx.resize(p->knox+8);
    cfy.resize(p->knoy+8);
    cfz.resize(p->knoz+8);
    
    isfx.resize(p->knox+8);
    isfy.resize(p->knoy+8);
    isfz.resize(p->knoz+8);
    
    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);
               
    iniflag=1;    
    }              
                      
}

void weno_nug_func::dsdiffx(slice &f, slice &dq)
{
    // faces i-3 .. i+2 of the cells 0 .. knox-1
    for(int ii=-3; ii<p->knox+2; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    dq(ii,jj) = (f(ii+1,jj)-f(ii,jj))/p->DXP[ii+marge];
}

void weno_nug_func::dsdiffy(slice &f, slice &dq)
{
    // faces j-3 .. j+2 of the cells 0 .. knoy-1
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-3; jj<p->knoy+2; ++jj)
    dq(ii,jj) = (f(ii,jj+1)-f(ii,jj))/p->DYP[jj+marge];
}

int weno_nug_func::iniflag(0);
