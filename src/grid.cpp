/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Alexander Hanke (@AlHanke), Hans Bihs
--------------------------------------------------------------------*/

#include "grid.h"
#include "lexer.h"
#include "ghostcell.h"
#include "resize.h"

void grid::assign_margin()
{
    imax=knox+2*margin;
    jmax=knoy+2*margin;
    kmax=knoz+2*margin;
    kmaxF=knoz+1+2*margin;

    imin=-margin;
    jmin=-margin;
    kmin=-margin;
}

/*!
    * @brief The sigma_coord_ini method initializes the sigma coordinates for vertical grids in REEF3D.
    * It calculates the normalized vertical coordinates (ZN) based on the physical vertical coordinates and the total depth of the grid.
    * The method also allocates memory for the ZSN and ZSP arrays, which store the sigma coordinates at nodal and cell-centered locations, respectively.
*/
void grid::sigma_coord_ini()
{
    double L, ZN0temp;

    L = ZN[gknoz+marge] - ZN[0+marge];

    ZN0temp = ZN[0+marge];

    for(k=0;k<max_k+1;++k)
    ZN[k] = (ZN[k]-ZN0temp)/L;

    resize_class resizer;
    resizer.Darray(ZSN,imax*jmax*(kmax+1));
    resizer.Darray(ZSP,imax*jmax*kmax);
}

/*!
    * @brief The gridspacing method calculates the grid spacing for the computational grid in REEF3D.
    * It initializes the cell-centered coordinates (XP, YP, ZP) and the grid spacing (DXN, DYN, DZN for nodal spacing and DXP, DYP, DZP for cell-centered spacing).
    * It also calculates the average grid spacing (DXM) and the minimum grid spacing in x and y directions (DXD, DYD).
    * The method uses a resize_class to allocate memory for the necessary arrays and performs global reductions to compute the average and minimum grid spacings across all MPI ranks.
*/
void grid::gridspacing(lexer* p, ghostcell *pgc)
{
    resize_class resizer;

    // Consolidate XN, YN, ZN into a single array for better cache performance
    increment::max_i = gknox + 2*marge;
    increment::max_j = gknoy + 2*marge;
    increment::max_k = gknoz + 2*marge;

    // x direction
    std::vector<int> transfer_number(amrex::ParallelDescriptor::NProcs(), 0);
    int local_size = knox;
    if(p->nb1==-2) local_size += marge;
    if(p->nb4==-2) local_size += marge+1;
    if(p->nb3!=-2 || p->nb5!=-2) local_size = 0;
    MPI_Allgather(&local_size, 1, MPI_INT, transfer_number.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    std::vector<int> offsets(amrex::ParallelDescriptor::NProcs(), 0);
    int local_offset = marge + p->origin_i;
    if(p->nb1==-2) local_offset -= marge;
    MPI_Allgather(&local_offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    std::vector<double> temp(max_i+1,0);
    int local_start = marge;
    if(p->nb1==-2) local_start -= marge;

    MPI_Allgatherv(XN.data()+local_start, local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    XN.clear(); XN.resize(max_i+1,0);
    std::copy(temp.begin(), temp.end(), XN.begin());

    MPI_Allgatherv(RN+local_start, local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    delete [] RN; RN = new double[max_i+1];
    std::copy(temp.begin(), temp.end(), RN);

    // y direction
    local_size = knoy;
    if(p->nb3==-2) local_size += marge;
    if(p->nb2==-2) local_size += marge+1;
    if(p->nb1!=-2 || p->nb5!=-2) local_size = 0;
    MPI_Allgather(&local_size, 1, MPI_INT, transfer_number.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    local_offset = marge + p->origin_j;
    if(p->nb3==-2) local_offset -= marge;
    MPI_Allgather(&local_offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    temp.clear();
    temp.resize(max_j+1,0);
    local_start = marge;
    if(p->nb3==-2) local_start -= marge;
    MPI_Allgatherv(YN.data(), local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    YN.clear(); YN.resize(max_j+1,0);
    std::copy(temp.begin(), temp.end(), YN.begin());

    MPI_Allgatherv(SN+local_start, local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    delete [] SN; SN = new double[max_j+1];
    std::copy(temp.begin(), temp.end(), SN);

    // z direction
    local_size = knoz;
    if(p->nb5==-2) local_size += marge;
    if(p->nb6==-2) local_size += marge+1;
    if(p->nb1!=-2 || p->nb3!=-2) local_size = 0;
    MPI_Allgather(&local_size, 1, MPI_INT, transfer_number.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    local_offset = marge + p->origin_k;
    if(p->nb5==-2) local_offset -= marge;
    MPI_Allgather(&local_offset, 1, MPI_INT, offsets.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());

    temp.clear();
    temp.resize(max_k+1,0);
    local_start = marge;
    if(p->nb5==-2) local_start -= marge;

    MPI_Allgatherv(ZN.data(), local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    ZN.clear(); ZN.resize(max_k+1,0);
    std::copy(temp.begin(), temp.end(), ZN.begin());

    MPI_Allgatherv(TN+local_start, local_size, MPI_DOUBLE, temp.data(), transfer_number.data(), offsets.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
    delete [] TN; TN = new double[max_k+1];
    std::copy(temp.begin(), temp.end(), TN);

    if(p->G2==1)
    sigma_coord_ini();

    // Derived coordinates and spacings
    XP.resize(max_i);
    YP.resize(max_j);
    ZP.resize(max_k);

    double *RP,*SP,*TP; // Temporary arrays
    resizer.Darray(RP,max_i);
    resizer.Darray(SP,max_j);
    resizer.Darray(TP,max_k);

    DXN.resize(max_i);
    DYN.resize(max_j);
    DZN.resize(max_k);

    DXP.resize(max_i);
    DYP.resize(max_j);
    DZP.resize(max_k);

    resizer.Darray(DRDXN,max_i);
    resizer.Darray(DSDYN,max_j);
    resizer.Darray(DTDZN,max_k);

    resizer.Darray(DRDXP,max_i);
    resizer.Darray(DSDYP,max_j);
    resizer.Darray(DTDZP,max_k);

    // XP,YP,ZP
    for(i=0;i<max_i;++i)
    XP[i] = 0.5*(XN[i]+XN[i+1]);

    for(j=0;j<max_j;++j)
    YP[j] = 0.5*(YN[j]+YN[j+1]);

    for(k=0;k<max_k;++k)
    ZP[k] = 0.5*(ZN[k]+ZN[k+1]);

    // RP,SP,TP
    for(i=0;i<max_i;++i)
    RP[i] = 0.5*(RN[i]+RN[i+1]);

    for(j=0;j<max_j;++j)
    SP[j] = 0.5*(SN[j]+SN[j+1]);

    for(k=0;k<max_k;++k)
    TP[k] = 0.5*(TN[k]+TN[k+1]);

    //dx
    for(i=0;i<max_i;++i)
    DXN[i] = XN[i+1]-XN[i];

    for(j=0;j<max_j;++j)
    DYN[j] = YN[j+1]-YN[j];

    for(k=0;k<max_k;++k)
    DZN[k] = ZN[k+1]-ZN[k];

    // dxn
    for(i=0;i<max_i;++i)
    DXP[i] = XP[i+1] - XP[i];

    for(j=0;j<max_j;++j)
    DYP[j] = YP[j+1] - YP[j];

    for(k=0;k<max_k;++k)
    DZP[k] = ZP[k+1] - ZP[k];

    // transformation
    // 1st derivative
    for(i=-1+marge;i<gknox+2+marge;++i)
    DRDXN[i] =  (-RN[i+2] + 8.0*RN[i+1] - 8.0*RN[i-1] + RN[i-2])
                /(-XN[i+2] + 8.0*XN[i+1] - 8.0*XN[i-1] + XN[i-2]);

    for(j=-1+marge;j<gknoy+2+marge;++j)
    DSDYN[j] =  (-SN[j+2] + 8.0*SN[j+1] - 8.0*SN[j-1] + SN[j-2])
                /(-YN[j+2] + 8.0*YN[j+1] - 8.0*YN[j-1] + YN[j-2]);

    for(k=-1+marge;k<gknoz+2+marge;++k)
    DTDZN[k] =  (-TN[k+2] + 8.0*TN[k+1] - 8.0*TN[k-1] + TN[k-2])
                /(-ZN[k+2] + 8.0*ZN[k+1] - 8.0*ZN[k-1] + ZN[k-2]);

    for(i=-1+marge;i<gknox+1+marge;++i)
    DRDXP[i] =  (-RP[i+2] + 8.0*RP[i+1] - 8.0*RP[i-1] + RP[i-2])
                /(-XP[i+2] + 8.0*XP[i+1] - 8.0*XP[i-1] + XP[i-2]);

    for(j=-1+marge;j<gknoy+1+marge;++j)
    DSDYP[j] =  (-SP[j+2] + 8.0*SP[j+1] - 8.0*SP[j-1] + SP[j-2])
                /(-YP[j+2] + 8.0*YP[j+1] - 8.0*YP[j-1] + YP[j-2]);

    for(k=-1+marge;k<gknoz+1+marge;++k)
    DTDZP[k] =  (-TP[k+2] + 8.0*TP[k+1] - 8.0*TP[k-1] + TP[k-2])
                /(-ZP[k+2] + 8.0*ZP[k+1] - 8.0*ZP[k-1] + ZP[k-2]);

    delete [] RN; RN = nullptr;
    delete [] SN; SN = nullptr;
    delete [] TN; TN = nullptr;
    delete [] RP; RP = nullptr;
    delete [] SP; SP = nullptr;
    delete [] TP; TP = nullptr;

    // Average grid spacing
    DXM = DXD = DYD = 0.0;

    int count=0;
    int xcount=0;
    int ycount=0;

    for(i=marge+origin_i;i<origin_i+marge+knox;++i)
    {
        DXM += DXP[i];
        DXD += DXP[i];
        ++count;
        ++xcount;
    }

    if(j_dir==1)
    for(j=marge+origin_j;j<origin_j+marge+knoy;++j)
    {
        DXM += DYP[j];
        DYD += DYP[j];
        ++count;
        ++ycount;
    }

    for(k=marge+origin_k;k<origin_k+marge+knoz;++k)
    {
        DXM += DZP[k];
        ++count;
    }

    count = pgc->globalisum(count);
    xcount = pgc->globalisum(xcount);
    ycount = pgc->globalisum(ycount);

    DXM = pgc->globalsum(DXM);
    DXD = pgc->globalsum(DXD);
    DYD = pgc->globalsum(DYD);

    DXM /= double(count);
    DXD /= double(xcount);
    DYD /= double(ycount);

    DXM = pgc->globalmin(DXM);
    DXD = pgc->globalmin(DXD);
    DYD = pgc->globalmin(DYD);
}
