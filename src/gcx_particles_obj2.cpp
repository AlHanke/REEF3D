/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include "ghostcell.h"
#include "lexer.h"
#include "particles_obj2.h"

#include <vector>

void ghostcell::para_particles_obj2(lexer* p, particles_obj2 s[6], particles_obj2 r[6])
{
    // Setup amount information
    int send[6], recv[6];
    for (int n=0;n<6;n++)
    {
        send[n]=s[n].particles.size();
        recv[n]=0;
    }

    // Transfer amount information
    if(p->nb1>=0)
    {
        MPI_Isend(&send[0],1,MPI_INT,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(&recv[0],1,MPI_INT,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->nb2>=0)
    {
        MPI_Isend(&send[1],1,MPI_INT,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(&recv[1],1,MPI_INT,p->nb2,tag3,mpi_comm,&rreq2);
    }

    if(p->nb3>=0)
    {
        MPI_Isend(&send[2],1,MPI_INT,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(&recv[2],1,MPI_INT,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->nb4>=0)
    {
        MPI_Isend(&send[3],1,MPI_INT,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(&recv[3],1,MPI_INT,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(&send[4],1,MPI_INT,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(&recv[4],1,MPI_INT,p->nb5,tag6,mpi_comm,&rreq5);
    }

    if(p->nb6>=0)
    {
        MPI_Isend(&send[5],1,MPI_INT,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(&recv[5],1,MPI_INT,p->nb6,tag5,mpi_comm,&rreq6);
    }

    gcwait(p);

    for (int n=0;n<6;n++)
        r[n].particles.resize(recv[n]);
    
    // Transfer data
    if(p->nb1>=0)
    {
        MPI_Isend(s[0].particles.data(),send[0],MPI_PARTICLE,p->nb1,tag1,mpi_comm,&sreq1);
        MPI_Irecv(r[0].particles.data(),recv[0],MPI_PARTICLE,p->nb1,tag4,mpi_comm,&rreq1);
	}

    if(p->nb2>=0)
    {
        MPI_Isend(s[1].particles.data(),send[1],MPI_PARTICLE,p->nb2,tag2,mpi_comm,&sreq2);
        MPI_Irecv(r[1].particles.data(),recv[1],MPI_PARTICLE,p->nb2,tag3,mpi_comm,&rreq2);
	}

    if(p->nb3>=0)
    {
        MPI_Isend(s[2].particles.data(),send[2],MPI_PARTICLE,p->nb3,tag3,mpi_comm,&sreq3);
        MPI_Irecv(r[2].particles.data(),recv[2],MPI_PARTICLE,p->nb3,tag2,mpi_comm,&rreq3);
	}

    if(p->nb4>=0)
    {
        MPI_Isend(s[3].particles.data(),send[3],MPI_PARTICLE,p->nb4,tag4,mpi_comm,&sreq4);
        MPI_Irecv(r[3].particles.data(),recv[3],MPI_PARTICLE,p->nb4,tag1,mpi_comm,&rreq4);
	}

    if(p->nb5>=0)
    {
        MPI_Isend(s[4].particles.data(),send[4],MPI_PARTICLE,p->nb5,tag5,mpi_comm,&sreq5);
        MPI_Irecv(r[4].particles.data(),recv[4],MPI_PARTICLE,p->nb5,tag6,mpi_comm,&rreq5);
	}

    if(p->nb6>=0)
    {
        MPI_Isend(s[5].particles.data(),send[5],MPI_PARTICLE,p->nb6,tag6,mpi_comm,&sreq6);
        MPI_Irecv(r[5].particles.data(),recv[5],MPI_PARTICLE,p->nb6,tag5,mpi_comm,&rreq6);
	}

    gcwait(p);
}