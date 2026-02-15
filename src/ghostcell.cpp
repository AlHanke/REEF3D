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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include<sstream>

ghostcell::ghostcell(int& argc, char **argv, lexer *p) : tag1(1),tag2(2),tag3(3),tag4(4),tag5(5),tag6(6)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_dup(MPI_COMM_WORLD,&mpi_comm);

    MPI_Comm_rank(mpi_comm,&p->mpirank);
    MPI_Comm_size(mpi_comm,&p->mpi_size);

    ghostcell::p=p;

    if(p->mpi_size==1)
        do_comms = false;
}

void ghostcell::gc_ini(lexer* p)
{
    margin=p->margin;
    paramargin=p->margin;
    gamma=p->B29;

    if(p->B23==1)
    orderdir=2;
    else if(p->B23==2)
    orderdir=3;

    if(p->B20==1)
    {
        gclabel_u=bc_labels::NEUMANN;
        gclabel_v=bc_labels::NEUMANN;
        gclabel_w=bc_labels::NEUMANN;
    }
    else if(p->B20==2)
    {
        gclabel_u=bc_labels::NOSLIP;
        gclabel_v=bc_labels::NOSLIP;
        gclabel_w=bc_labels::NOSLIP;
    }
    else if(p->B20==3)
    {
        gclabel_u=bc_labels::DIRICHLET_ORTH;
        gclabel_v=bc_labels::DIRICHLET_ORTH;
        gclabel_w=bc_labels::DIRICHLET_ORTH;
    }
    else if(p->B20==4)
    {
        gclabel_u=bc_labels::NOSLIP;
        gclabel_v=bc_labels::NOSLIP;
        gclabel_w=bc_labels::NOSLIP;
    }

    gclabel_u_orth=bc_labels::DIRICHLET_ORTH;
    gclabel_v_orth=bc_labels::DIRICHLET_ORTH;
    gclabel_w_orth=bc_labels::DIRICHLET_ORTH;

    // for reflective BC
    if(p->B23==2)
    {
        gclabel_u=bc_labels::DIRICHLET_PARA_REFLECT;
        gclabel_v=bc_labels::DIRICHLET_PARA_REFLECT;
        gclabel_w=bc_labels::DIRICHLET_PARA_REFLECT;

        gclabel_u_orth=bc_labels::DIRICHLET_ORTH_REFLECT;
        gclabel_v_orth=bc_labels::DIRICHLET_ORTH_REFLECT;
        gclabel_w_orth=bc_labels::DIRICHLET_ORTH_REFLECT;
    }

    awa_label=false;
    if(p->B99>=3)
    awa_label=true;


    if(p->B75==1)
    {
        gclabel_u_out=bc_labels::NEUMANN;
        gclabel_v_out=bc_labels::NEUMANN;
        gclabel_w_out=bc_labels::NEUMANN;
    }
    else if(p->B75==2)
    {
        gclabel_u_out=bc_labels::OUTFLOW;
        gclabel_v_out=bc_labels::OUTFLOW;
        gclabel_w_out=bc_labels::OUTFLOW;
    }
    else if(p->B75==3)
    {
        gclabel_u_out=bc_labels::NONE;
        gclabel_v_out=bc_labels::OUTFLOW;
        gclabel_w_out=bc_labels::OUTFLOW;
    }

    gclabel_outflow=true;
    if(p->B60==3||p->B60==4)
    gclabel_outflow=false;

    gclabel_u_in=bc_labels::DIRICHLET_ORTH;
    gclabel_v_in=bc_labels::DIRICHLET_ORTH;
    gclabel_w_in=bc_labels::DIRICHLET_ORTH;
    gclabel_lsm_in=bc_labels::NEUMANN;

    if(p->I230>0 || p->B98>=3 || p->B60>0)
    {
        gclabel_u_in=bc_labels::NONE;
        gclabel_v_in=bc_labels::NONE;
        gclabel_w_in=bc_labels::NONE;
        gclabel_lsm_in=bc_labels::NONE;
    }

    // pressure bc labels
    gclabel_press_in=bc_labels::NEUMANN;
    if(p->B76==2 || p->B76==3)
    gclabel_press_in=bc_labels::NONE;

    // pressure inflow
    pressin_label=false;
    if(p->B76!=1)
    pressin_label=true;

    // pressure outflow
    pressout_label=false;
    if(p->B77==1 || p->B77==10)
    pressout_label=true;

    // sflow slip/no-slip
    if(p->A217==1 && p->A10==2)
    {
        gclabel_u=bc_labels::NEUMANN;
        gclabel_v=bc_labels::NEUMANN;
    }
    else if(p->A217==2 && p->A10==2)
    {
        gclabel_u=bc_labels::NOSLIP;
        gclabel_v=bc_labels::NOSLIP;
    }

    for(m=0;m<15;m++)
    {
        y[m]=0.0;
    }
}

void ghostcell::fdm2D_update(fdm2D *bb)
{
    b=bb;
}

void ghostcell::fdm_fnpf_update(fdm_fnpf *cc)
{
    c=cc;
}

void ghostcell::fdm_nhf_update(fdm_nhf *dd)
{
    d=dd;
}

void ghostcell::fdm_update(fdm *aa)
{
    a=aa;
}

