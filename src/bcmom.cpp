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

#include "bcmom.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "gcb_sl_list.h"

bcmom::bcmom(lexer* p): surftens(p), roughness()
{
}

void bcmom::bcmom_start(fdm* a, lexer* p, ghostcell *pgc, turbulence *pturb, field& b, int gcval)
{
    wall_laws(p,a,b,gcval);
    surftens::surface_tension(p,a,a->phi,gcval);
}

void bcmom::wall_laws(lexer* p, fdm* a, field& b, int gcval)
{
    if(p->B10!=0)
    {
        int q;

        if(gcval==10)
        {
            QGC1LOOP
            {
                auto &gcb_entry = p->gcb1[p->level][q];
                if(gcb_entry.bc==21 && gcb_entry.cs!=X_NEG && gcb_entry.cs!=X_POS)
                    wall_law_u(p,a,b,gcb_entry,p->level);
            }

            QGCDF1LOOP
                wall_law_u(p,a,b,p->gcdf1[p->level][q],p->level);
        }
        else if(gcval==11 && p->j_dir==1)
        {
            QGC2LOOP
            {
                auto &gcb_entry = p->gcb2[p->level][q];
                if(gcb_entry.bc==21 && gcb_entry.cs!=Y_POS && gcb_entry.cs!=Y_NEG)
                    wall_law_v(p,a,b,gcb_entry,p->level);
            }

            QGCDF2LOOP
                wall_law_v(p,a,b,p->gcdf2[p->level][q],p->level);
        }
        else if(gcval==12)
        {
            QGC3LOOP
            {
                auto &gcb_entry = p->gcb3[p->level][q];
                if(gcb_entry.bc==21 && gcb_entry.cs!=Z_NEG && gcb_entry.cs!=Z_POS)
                    wall_law_w(p,a,b,gcb_entry,p->level);
            }

            QGCDF3LOOP
                wall_law_w(p,a,b,p->gcdf3[p->level][q],p->level);

        }
    }
}

template<typename gcb_entry_t>
void bcmom::wall_law_u(lexer* p, fdm* a, field& b, gcb_entry_t &gcb_entry, int lev)
{
    int i = gcb_entry.i;
    int j = gcb_entry.j;
    int k = gcb_entry.k;

    GCB_TILE(gcb_entry, lev);

    int cs = gcb_entry.cs;

    if(cs==Y_POS || cs==Y_NEG)
        deltaZ = p->DYN[JP];
    else if(cs==Z_NEG || cs==Z_POS)
        deltaZ = p->DZN[KP];

    z0 = 0.5*deltaZ;

    ks = ks_val(p,a,i,j,k,cs);

    if(30.0*z0<ks)
        z0 = ks/30.0;

    uplus = (1.0/kappa)*log(30.0*(z0/ks));

    a->F(i,j,k) -= ((fabs(a->u(i,j,k))*a->u(i,j,k))/(uplus*uplus*deltaZ));
}

template<typename gcb_entry_t>
void bcmom::wall_law_v(lexer* p, fdm* a, field& b, gcb_entry_t &gcb_entry, int lev)
{
    int i = gcb_entry.i;
    int j = gcb_entry.j;
    int k = gcb_entry.k;

    GCB_TILE(gcb_entry, lev);

    int cs = gcb_entry.cs;

    if(cs==X_NEG || cs==X_POS)
        deltaZ = p->DXN[IP];
    else if(cs==Z_NEG || cs==Z_POS)
        deltaZ = p->DZN[KP];

    z0 = 0.5*deltaZ;

    ks = ks_val(p,a,i,j,k,cs);

    if(30.0*z0<ks)
        z0=ks/30.0;

    uplus = (1.0/kappa)*log(30.0*(z0/ks));

    a->G(i,j,k) -= ((fabs(a->v(i,j,k))*a->v(i,j,k))/(uplus*uplus*deltaZ));
}

template<typename gcb_entry_t>
void bcmom::wall_law_w(lexer* p, fdm* a, field& b, gcb_entry_t &gcb_entry, int lev)
{
    int i = gcb_entry.i;
    int j = gcb_entry.j;
    int k = gcb_entry.k;

    GCB_TILE(gcb_entry, lev);

    int cs = gcb_entry.cs;

    if(cs==X_NEG || cs==X_POS)
        deltaZ = p->DXN[IP];
    else if(cs==Y_POS || cs==Y_NEG)
        deltaZ = p->DYN[JP];

    z0 = 0.5*deltaZ;

    ks = ks_val(p,a,i,j,k,cs);

    if(30.0*z0<ks)
        z0 = ks/30.0;

    uplus = (1.0/kappa)*log(30.0*(z0/ks));

    a->H(i,j,k) -= ((fabs(a->w(i,j,k))*a->w(i,j,k))/(uplus*uplus*deltaZ));
}
