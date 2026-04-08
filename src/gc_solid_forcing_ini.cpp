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
Authors: Hans Bihs, Tobias Martin, Ahmet Soydan
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"heaviside_ls.h"

void ghostcell::solid_forcing_ini(lexer *p, fdm *a)
{
    // ghostcell update
    gcdf_update(p,a);

    // Initialise floating fields
    ULOOP
    a->fbh1(i,j,k) = Hsolidface(p,a->solid,a->topo,1,0,0);

    VLOOP
    a->fbh2(i,j,k) = Hsolidface(p,a->solid,a->topo,0,1,0);

    WLOOP
    a->fbh3(i,j,k) = Hsolidface(p,a->solid,a->topo,0,0,1);

    LOOP
    a->fbh4(i,j,k) = Hsolidface(p,a->solid,a->topo,0,0,0);

    start1(p,a->fbh1,10);
    start2(p,a->fbh2,11);
    start3(p,a->fbh3,12);
    start4(p,a->fbh4,40);

    double dirac;
    LOOP
    {
        double psi;
        if(!p->j_dir)
        psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
        else
        psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        if(fabs(MIN(a->solid(i,j,k),a->topo(i,j,k)))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(MIN(a->solid(i,j,k),a->topo(i,j,k))))/psi));
        else
        dirac = 0.0;

        a->fbh5(i,j,k) =  1.0-MIN(dirac,1.0);
    }
}

template<typename GenericFieldConst>
inline double ghostcell::Hsolidface(lexer *p, GenericFieldConst& solid, GenericFieldConst& topo, int aa, int bb, int cc)
{
    if(p->topoforcing==0 && p->solidread==0)
    return 0.0;

    double psi;
    if(!p->j_dir)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
    else
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    double phival_sf;
    if(p->topoforcing>0 && p->solidread==0)
    phival_sf = 0.5*(topo(i,j,k) + topo(i+aa,j+bb,k+cc));
    else if(p->topoforcing==0 && p->solidread>0)
    phival_sf = 0.5*(solid(i,j,k) + solid(i+aa,j+bb,k+cc));
    else if(p->topoforcing>0 && p->solidread>0)
    phival_sf = std::min(0.5*(solid(i,j,k) + solid(i+aa,j+bb,k+cc)), 0.5*(topo(i,j,k) + topo(i+aa,j+bb,k+cc)));
    // else == if(p->topoforcing==0 && p->solidread==0) is covered in the beginning of the function

    return heaviside_ls(-phival_sf, psi);
}
