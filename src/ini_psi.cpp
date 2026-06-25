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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::inipsi(lexer* p, fdm *a, ghostcell* pgc)
{
    double psim;
    int count;
    
    if(p->j_dir==0)        
    p->psi = p->F45*(1.0/2.0)*(p->DXM+p->DZM);
        
    if(p->j_dir==1)
    p->psi = p->F45*(1.0/3.0)*(p->DXM+p->DYM+p->DZM);
    
    
    if(p->B90>0 || p->B60>0)
    {
    // psi
        count=0;
        psim=0.0;
        LOOP
        if(fabs(a->phi(i,j,k))<5.0*p->DZM)
        {
        psim += p->DZN[KP];
        ++count;
        }
        
        count=pgc->globalisum(count);
        psim=pgc->globalsum(psim);
        
        p->psi = p->F45*psim/double(count);

    }

    // Per-level interface half-width: the smoothed-Heaviside smearing width must scale with
    // the local cell size, else a finer AMR level over-smears the interface (and the coarse
    // roface no longer matches the average of the covering fine-cell roface values, breaking
    // density consistency across the C-F interface). psi is set from level-0 mean spacing, so
    // scale it down per level by that level's cell-size ratio. Level 0 (and single-level runs)
    // keep the exact original value.
    #if USE_AMREX
    double ratio;
    if(p->j_dir==0)
        ratio=0.5*(1.0/double(p->ref_vec[0]) + 1.0/double(p->ref_vec[2]))/1.0;
    else
        ratio=(1.0/double(p->ref_vec[0]) + 1.0/double(p->ref_vec[1]) + 1.0/double(p->ref_vec[2]))/3.0;
    double f=1.0; for(int n=0;n<p->nlevs-1;++n) f*=ratio;
    p->psi *= f; // scale down the level-0 psi to match the finest level's cell size

    p->psi *= 0.5; // temp
    #endif
}