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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include <HYPRE_utilities.h>
#include<sstream>

ghostcell::ghostcell(int& argc, char **argv, lexer *p) : margin(p->margin)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_dup(MPI_COMM_WORLD,&mpi_comm);

    MPI_Comm_rank(mpi_comm,&p->mpirank);
    MPI_Comm_size(mpi_comm,&p->mpi_size);

    HYPRE_Initialize();

    ghostcell::p=p;

    if(p->mpi_size==1)
    do_comms = false;
}

void ghostcell::gc_ini(lexer* p)
{
    paramargin=margin;
    gamma=p->B29;

    if(p->B23==1)
        orderdir=2;
    else if(p->B23==2)
        orderdir=3;

    if(p->B20==1)
    {
        gclabel_u=4;
        gclabel_v=4;
        gclabel_w=4;

        gclabel_utopo=4;
        gclabel_vtopo=4;
        gclabel_wtopo=4;
    }
    else if(p->B20==2)
    {
        gclabel_u=5;
        gclabel_v=5;
        gclabel_w=5;

        gclabel_utopo=5;
        gclabel_vtopo=5;
        gclabel_wtopo=5;
    }
    else if(p->B20==3)
    {
        gclabel_u=2;
        gclabel_v=2;
        gclabel_w=2;

        gclabel_utopo=2;
        gclabel_vtopo=2;
        gclabel_wtopo=2;
    }
    else if(p->B20==4)
    {
        gclabel_u=5;
        gclabel_v=5;
        gclabel_w=5;

        gclabel_utopo=4;
        gclabel_vtopo=4;
        gclabel_wtopo=4;
    }

    gclabel_u_orth=1;
    gclabel_v_orth=1;
    gclabel_w_orth=1;

    // for reflective BC
    if(p->B23==2)
    {
        gclabel_u=12;
        gclabel_v=12;
        gclabel_w=12;

        gclabel_u_orth=11;
        gclabel_v_orth=11;
        gclabel_w_orth=11;
    }

    awa_label=false;
    if(p->B99>=3)
        awa_label=true;


    if(p->B75==1)
    {
        gclabel_u_out=4;
        gclabel_v_out=4;
        gclabel_w_out=4;
    }
    else if(p->B75==2)
    {
        gclabel_u_out=6;
        gclabel_v_out=6;
        gclabel_w_out=6;
    }
    else if(p->B75==3)
    {
        gclabel_u_out=0;
        gclabel_v_out=6;
        gclabel_w_out=6;
    }

    gclabel_outflow=true;
    if(p->B60==3||p->B60==4)
        gclabel_outflow=false;

    gclabel_u_in=1;
    gclabel_v_in=1;
    gclabel_w_in=1;
    gclabel_lsm_in=4;

    if(p->I230>0 || p->B98>=3 || p->B60>0)
    {
        gclabel_u_in=0;
        gclabel_v_in=0;
        gclabel_w_in=0;
        gclabel_lsm_in=0;
    }

    // pressure bc labels
    gclabel_press_in=4;
    if(p->B76==2 || p->B76==3)
        gclabel_press_in=0;

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
        gclabel_u=4;
        gclabel_v=4;
    }
    else if(p->A217==2 && p->A10==2)
    {
        gclabel_u=5;
        gclabel_v=5;
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

void ghostcell::final(bool error)
{
    HYPRE_Finalize();
    if(cart_comm != MPI_COMM_NULL)
        MPI_Comm_free(&cart_comm);
    if(mpi_comm != MPI_COMM_NULL)
        MPI_Comm_free(&mpi_comm);
    MPI_Finalize();
    exit(error);
}
