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

#include"ghostcell.h"
#include"lexer.h"
#include"time.h"

double ghostcell::globalsum(double sendsum)
{
    MPI_Allreduce(&sendsum,&recvsum,1,MPI_DOUBLE,MPI_SUM,mpi_comm);
    return recvsum;
}

int ghostcell::globalisum(int sendisum)
{
    MPI_Allreduce(&sendisum,&recvisum,1,MPI_INT,MPI_SUM,mpi_comm);
    return recvisum;
}

double ghostcell::globalmin(double sendmin)
{
    MPI_Allreduce(&sendmin,&recvmin,1,MPI_DOUBLE,MPI_MIN,mpi_comm);
    return recvmin;
}

double ghostcell::globalmax(double sendmax)
{
    MPI_Allreduce(&sendmax,&recvmax,1,MPI_DOUBLE,MPI_MAX,mpi_comm);
    return recvmax;
}

int ghostcell::globalimin(int sendimin)
{
    MPI_Allreduce(&sendimin,&recvimin,1,MPI_INT,MPI_MIN,mpi_comm);
    return recvimin;
}

int ghostcell::globalimax(int sendimax)
{
    MPI_Allreduce(&sendimax,&recvimax,1,MPI_INT,MPI_MAX,mpi_comm);
    return recvimax;
}

double ghostcell::timesync(double t)
{
    MPI_Bcast(&t,1,MPI_DOUBLE,0,mpi_comm);
    return t;
}

void ghostcell::globalctrl(lexer* p)
{
    MPI_Bcast(p->ictrl,p->ctrlsize,MPI_INT,0,mpi_comm);
    MPI_Bcast(p->dctrl,p->ctrlsize,MPI_DOUBLE,0,mpi_comm);
}

double ghostcell::timer()
{
    double t=0.0;
    t=MPI_Wtime();
    return t;
}


int ghostcell::Bcast(void *buffer, int count, MPI_Datatype datatype)
{
    return MPI_Bcast(buffer,count,datatype,0,mpi_comm);;
}
int ghostcell::File_open_createWriteOnly(MPI_File *fh, const char *filename)
{
    return MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, fh);
}
int ghostcell::File_set_size(MPI_File fh, MPI_Offset size)
{
    return MPI_File_set_size(fh, size);
}
int ghostcell::File_write_at_char(MPI_File fh, MPI_Offset offset, const void *buf, int count)
{
    return MPI_File_write_at(fh, offset, buf, count, MPI_CHAR, MPI_STATUS_IGNORE);
}
int ghostcell::File_write_at_all_float(MPI_File fh, MPI_Offset offset, const void *buf, int count)
{
    return MPI_File_write_at_all(fh, offset, buf, count, MPI_FLOAT, MPI_STATUS_IGNORE);
}
int ghostcell::File_write_at_all_char(MPI_File fh, MPI_Offset offset, const void *buf, int count)
{
    return MPI_File_write_at_all(fh, offset, buf, count, MPI_CHAR, MPI_STATUS_IGNORE);
}
int ghostcell::File_close(MPI_File *fh)
{
    return MPI_File_close(fh);
}