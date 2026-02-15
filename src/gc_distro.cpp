/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Alexander Hanke

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "ghostcell.h"
#include "field.h"

void ghostcell::gcdistro(field& f, int ii, int jj, int kk, double dist, bc_labels bc_label, int cs)
{
    i=ii;
    j=jj;
    k=kk;

    switch(bc_label)
    {
        case bc_labels::DIRICHLET_ORTH:
            dirichlet_ortho(f,dist,cs);
            break;
        case bc_labels::NEUMANN:
            neumann(f,cs);
            break;
        case bc_labels::NOSLIP:
            noslip(f,cs);
            break;
        case bc_labels::OUTFLOW:
            outflow(f,cs);
            break;
        case bc_labels::POTENTIAL:
            potentialbc(f,cs);
            break;
        case bc_labels::DIRICHLET_ORTH_REFLECT:
            dirichlet_ortho_reflect(f,dist,cs);
            break;
        case bc_labels::DIRICHLET_PARA_REFLECT:
            dirichlet_para_reflect(f,dist,cs);
            break;
        case bc_labels::HEATBC:
            heatbc(f,cs);
            break;
        default:
            break;
    }
}

void ghostcell::gcsldistro(lexer *p, slice &f, int ii, int jj, bc_labels bc_label, int cs)
{
    i=ii;
    j=jj;

    switch (bc_label)
    {
    case bc_labels::NEUMANN:
        gcsl_neumann(f,cs);
        break;
    case bc_labels::NEUMANN_X:
        gcsl_neumann_x(f,cs);
        break;
    case bc_labels::NEUMANN_HX:
        gcsl_neumann_hx(f,cs);
        break;
    case bc_labels::NEUMANN_HY:
        gcsl_neumann_hy(f,cs);
        break;
    case bc_labels::NOSLIP:
        gcsl_noslip(f,cs);
        break;
    case bc_labels::OUTFLOW:
        gcsl_outflow(f,cs);
        break;
    case bc_labels::SOMMERFELD:
        gcsl_sommerfeld(p,f,cs);
        break;
    case bc_labels::POTENTIAL:
        gcsl_potentialbc(p,f,cs);
        break;
    default:
        break;
    }
}