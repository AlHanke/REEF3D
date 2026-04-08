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

#include"lexer.h"
#include<algorithm>
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void lexer::read_grid()
{
    int i,n;
    int isurf,jsurf,ksurf,surfside,surfgroup,side1,side2,paraconum;
    int para_active;
    char name[100];
    int DM_M10;
    int iin;
    double ddn;

    int gcwall_count=0;
    gcin_count=0;
    gcout_count=0;
    gcpara1_count=0;
    gcpara2_count=0;
    gcpara3_count=0;
    gcpara4_count=0;
    gcpara5_count=0;
    gcpara6_count=0;
    gcparaco1_count=0;
    gcparaco2_count=0;
    gcparaco3_count=0;
    gcparaco4_count=0;
    gcparaco5_count=0;
    gcparaco6_count=0;

    const int padding = 6;
    sprintf(name,"DIVEMesh_Grid/grid-%0*i.dat",padding,mpirank+1);

    // open file------------
    ifstream grid(name, ios_base::binary);

    if(!grid)
    {
        if(mpirank==0)
        {
            cout<<endl;
            cout<<"!!! Could not open DIVEMesh grid file: "<<name<<" !"<<endl;
            cout<<"!!! please check the manual!"<<endl<<endl<<endl<<endl;
        }
        exit(1);
    }

    // read grid file-------------

    grid.read((char*)&iin, sizeof (int));
    DM_M10=iin;

    if(mpirank==0)
    if(DM_M10!=M10 || M10!=mpi_size || DM_M10!=mpi_size)
    {
        cout<<endl;
        cout<<"!!! Inconsistent M 10 parameter, needs to be the same in REEF3D and DIVEMesh !"<<endl;
        cout<<"mpi_size: "<<mpi_size<<" REEFD M10: "<<M10<<" DIVEMesh M10: "<<DM_M10<<endl;
        cout<<"!!! please check the manual!"<<endl<<endl<<endl<<endl;

        exit(1);
    }

    grid.read((char*)&iin, sizeof (int));
    knox=iin;
    grid.read((char*)&iin, sizeof (int));
    knoy=iin;
    grid.read((char*)&iin, sizeof (int));
    knoz=iin;

    grid.read((char*)&ddn, sizeof (double));
    dx=ddn;

    grid.read((char*)&ddn, sizeof (double));
    DXM=ddn;
    grid.read((char*)&ddn, sizeof (double));
    DYM=ddn;
    grid.read((char*)&ddn, sizeof (double));
    DZM=ddn;

    grid.read((char*)&ddn, sizeof (double));
    originx=ddn;
    grid.read((char*)&ddn, sizeof (double));
    originy=ddn;
    grid.read((char*)&ddn, sizeof (double));
    originz=ddn;

    grid.read((char*)&ddn, sizeof (double));
    endx=ddn;
    grid.read((char*)&ddn, sizeof (double));
    endy=ddn;
    grid.read((char*)&ddn, sizeof (double));
    endz=ddn;

    grid.read((char*)&ddn, sizeof (double));
    global_xmin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_ymin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_zmin=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_xmax=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_ymax=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_zmax=ddn;

    grid.read((char*)&iin, sizeof (int));
    gknox=iin;
    grid.read((char*)&iin, sizeof (int));
    gknoy=iin;
    grid.read((char*)&iin, sizeof (int));
    gknoz=iin;

    grid.read((char*)&iin, sizeof (int));
    origin_i=iin;
    grid.read((char*)&iin, sizeof (int));
    origin_j=iin;
    grid.read((char*)&iin, sizeof (int));
    origin_k=iin;

    grid.read((char*)&iin, sizeof (int));
    gcwall_count=iin;

    grid.read((char*)&iin, sizeof (int));
    gcpara1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara4_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara5_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcpara6_count=iin;

    grid.read((char*)&iin, sizeof (int));
    gcparaco1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco4_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco5_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcparaco6_count=iin;

    grid.read((char*)&iin, sizeof (int));
    gcslpara1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslpara4_count=iin;

    grid.read((char*)&iin, sizeof (int));
    gcslparaco1_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco2_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco3_count=iin;
    grid.read((char*)&iin, sizeof (int));
    gcslparaco4_count=iin;

    grid.read((char*)&iin, sizeof (int));
    nb1=iin;
    grid.read((char*)&iin, sizeof (int));
    nb2=iin;
    grid.read((char*)&iin, sizeof (int));
    nb3=iin;
    grid.read((char*)&iin, sizeof (int));
    nb4=iin;
    grid.read((char*)&iin, sizeof (int));
    nb5=iin;
    grid.read((char*)&iin, sizeof (int));
    nb6=iin;

    grid.read((char*)&iin, sizeof (int));
    mx=iin;
    grid.read((char*)&iin, sizeof (int));
    my=iin;
    grid.read((char*)&iin, sizeof (int));
    mz=iin;

    grid.read((char*)&iin, sizeof (int));
    bcside1=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside2=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside3=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside4=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside5=iin;
    grid.read((char*)&iin, sizeof (int));
    bcside6=iin;

    grid.read((char*)&iin, sizeof (int));
    periodic1=iin;
    grid.read((char*)&iin, sizeof (int));
    periodic2=iin;
    grid.read((char*)&iin, sizeof (int));
    periodic3=iin;

    grid.read((char*)&iin, sizeof (int));
    periodicX1=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX2=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX3=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX4=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX5=iin;
    grid.read((char*)&iin, sizeof (int));
    periodicX6=iin;

    grid.read((char*)&iin, sizeof (int));
    i_dir=iin;
    grid.read((char*)&iin, sizeof (int));
    j_dir=iin;
    grid.read((char*)&iin, sizeof (int));
    k_dir=iin;

    grid.read((char*)&iin, sizeof (int));
    P150=iin;

    grid.read((char*)&iin, sizeof (int));
    solidread=iin; // solid
    grid.read((char*)&iin, sizeof (int));
    toporead=iin; // topo
    grid.read((char*)&iin, sizeof (int));
    solid_gcb_est=iin;
    grid.read((char*)&iin, sizeof (int));
    topo_gcb_est=iin;

    grid.read((char*)&iin, sizeof (int));
    solid_gcbextra_est=iin;
    grid.read((char*)&iin, sizeof (int));
    topo_gcbextra_est=iin;
    grid.read((char*)&iin, sizeof (int));
    tot_gcbextra_est=iin;

    grid.read((char*)&iin, sizeof (int));
    cms_flag=iin;
 
    grid.read((char*)&ddn, sizeof (double));
    global_orig_x=ddn;
    grid.read((char*)&ddn, sizeof (double));
    global_orig_y=ddn;
    grid.read((char*)&ddn, sizeof (double));
    alpha_grid=ddn;

    // Refinement
    grid.read((char*)&iin, sizeof (int));
    #if USE_AMREX
    nlevs = iin;
    #endif

    // ---------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------

    topo_gcb_est*=4;

    gcb1_count=gcb2_count=gcb3_count=gcb4_count=gcwall_count;

    gcpara_sum=gcpara1_count+gcpara2_count+gcpara3_count+gcpara4_count+gcpara5_count+gcpara6_count;
    gcparaco1_count+gcparaco2_count+gcparaco3_count+gcparaco4_count+gcparaco5_count+gcparaco6_count;

    grid::assign_margin();

    flag4_grid = std::unique_ptr<int[]>(new int[imax*jmax*kmax]);
    std::fill_n(flag4_grid.get(), imax*jmax*kmax, -1);

    //if(solidread==1)
    flag_solid = std::make_unique<int[]>(imax*jmax*kmax);

    //if(toporead==1)
    flag_topo = std::make_unique<int[]>(imax*jmax*kmax);

    solidbed = std::make_unique<double[]>(imax*jmax);
    topobed = std::make_unique<double[]>(imax*jmax);
    Darray(bed,imax*jmax);
    Iarray(wet,imax*jmax);
    Iarray(deep,imax*jmax);
    Darray(WL,imax*jmax);
    data = std::make_unique<double[]>(imax*jmax*kmax);
    Iarray(flagslice1,imax*jmax);
    Iarray(flagslice2,imax*jmax);
    Iarray(flagslice4,imax*jmax);

    for(i=0;i<imax*jmax;++i)
    {
        flagslice1[i]=-10;
        flagslice2[i]=-10;
        flagslice4[i]=-10;
    }

    if(gcb4_count>0)
    {
        Iarray(gcb1, gcb1_count,6);
        Iarray(gcb2, gcb2_count,6);
        Iarray(gcb3, gcb3_count,6);
        Iarray(gcb4, gcb4_count,6);
    }

    gcpara1.resize(gcpara1_count);
    gcpara2.resize(gcpara2_count);
    gcpara3.resize(gcpara3_count);
    gcpara4.resize(gcpara4_count);
    gcpara5.resize(gcpara5_count);
    gcpara6.resize(gcpara6_count);

    gcparaco1.resize(gcparaco1_count);
    gcparaco2.resize(gcparaco2_count);
    gcparaco3.resize(gcparaco3_count);
    gcparaco4.resize(gcparaco4_count);
    gcparaco5.resize(gcparaco5_count);
    gcparaco6.resize(gcparaco6_count);

    // Slice allocation
    gcbsl1_count=gcbsl2_count=gcbsl4_count=1;

    Iarray(gcbsl1, gcbsl1_count,5);
    Iarray(gcbsl2, gcbsl2_count,5);
    Iarray(gcbsl4, gcbsl4_count,5);

    gcslpara1.resize(gcslpara1_count);
    gcslpara2.resize(gcslpara2_count);
    gcslpara3.resize(gcslpara3_count);
    gcslpara4.resize(gcslpara4_count);

    gcslparaco1.resize(gcslparaco1_count);
    gcslparaco2.resize(gcslparaco2_count);
    gcslparaco3.resize(gcslparaco3_count);
    gcslparaco4.resize(gcslparaco4_count);

    XN.resize(knox+2*marge+1,0);
    YN.resize(knoy+2*marge+1,0);
    ZN.resize(knoz+2*marge+1,0);


    // ---------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------

    //  Flag
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
        grid.read((char*)&iin, sizeof (int));
        flag4_grid[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=iin;
    }

    // Nodes XYZ
    for(i=-marge;i<knox+marge+1;++i)
    {
        grid.read((char*)&ddn, sizeof (double));
        XN[i + marge]=ddn;
    }

    for(j=-marge;j<knoy+1+marge;++j)
    {
        grid.read((char*)&ddn, sizeof (double));
        YN[j + marge]=ddn;
    }

    for(k=-marge;k<knoz+1+marge;++k)
    {
        grid.read((char*)&ddn, sizeof (double));
        ZN[k + marge]=ddn;
    }

    //  Solid
    if(solidread==1)
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
        grid.read((char*)&ddn, sizeof (double));
        flag_solid[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=ddn;
    }

    //  Topo
    if(toporead==1)
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    for(k=0; k<knoz; ++k)
    {
        grid.read((char*)&ddn, sizeof (double));
        flag_topo[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=ddn;
    }

    //  GC Surfaces
    gcin_count=0;
    gcout_count=0;
    for(i=0; i<gcb4_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        grid.read((char*)&iin, sizeof (int));
        surfside=iin;

        grid.read((char*)&iin, sizeof (int));
        surfgroup=iin;

        gcb4[i][0]=isurf;
        gcb4[i][1]=jsurf;
        gcb4[i][2]=ksurf;
        gcb4[i][3]=surfside;
        gcb4[i][4]=surfgroup;

        if(surfgroup==1 || surfgroup==6)
            ++gcin_count;
        else if(surfgroup==2 || surfgroup==7)
            ++gcout_count;
    }

    Iarray(gcin, gcin_count,4);
    Iarray(gcout, gcout_count,4);

    //  Para Surfaces
    for(i=0; i<gcpara1_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara1[i][0]=isurf;
        gcpara1[i][1]=jsurf;
        gcpara1[i][2]=ksurf;
        gcpara1[i][3]=1;
    }

    for(i=0; i<gcpara2_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara2[i][0]=isurf;
        gcpara2[i][1]=jsurf;
        gcpara2[i][2]=ksurf;
        gcpara2[i][3]=1;
    }

    for(i=0; i<gcpara3_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara3[i][0]=isurf;
        gcpara3[i][1]=jsurf;
        gcpara3[i][2]=ksurf;
        gcpara3[i][3]=1;
    }

    for(i=0; i<gcpara4_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara4[i][0]=isurf;
        gcpara4[i][1]=jsurf;
        gcpara4[i][2]=ksurf;
        gcpara4[i][3]=1;
    }

    for(i=0; i<gcpara5_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara5[i][0]=isurf;
        gcpara5[i][1]=jsurf;
        gcpara5[i][2]=ksurf;
        gcpara5[i][3]=1;
    }

    for(i=0; i<gcpara6_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcpara6[i][0]=isurf;
        gcpara6[i][1]=jsurf;
        gcpara6[i][2]=ksurf;
        gcpara6[i][3]=1;
    }

    //  Para Corners
    for(i=0; i<gcparaco1_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco1[i][0]=isurf;
        gcparaco1[i][1]=jsurf;
        gcparaco1[i][2]=ksurf;
    }

    for(i=0; i<gcparaco2_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco2[i][0]=isurf;
        gcparaco2[i][1]=jsurf;
        gcparaco2[i][2]=ksurf;
    }

    for(i=0; i<gcparaco3_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco3[i][0]=isurf;
        gcparaco3[i][1]=jsurf;
        gcparaco3[i][2]=ksurf;
    }

    for(i=0; i<gcparaco4_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco4[i][0]=isurf;
        gcparaco4[i][1]=jsurf;
        gcparaco4[i][2]=ksurf;
    }

    for(i=0; i<gcparaco5_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco5[i][0]=isurf;
        gcparaco5[i][1]=jsurf;
        gcparaco5[i][2]=ksurf;
    }

    for(i=0; i<gcparaco6_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        grid.read((char*)&iin, sizeof (int));
        ksurf=iin;

        gcparaco6[i][0]=isurf;
        gcparaco6[i][1]=jsurf;
        gcparaco6[i][2]=ksurf;
    }

    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        grid.read((char*)&iin, sizeof (int));
        flagslice4[(i-imin)*jmax + (j-jmin)]=iin;
    }

    //  Paraslice Surfaces
    for(i=0; i<gcslpara1_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslpara1[i][0]=isurf;
        gcslpara1[i][1]=jsurf;
    }

    for(i=0; i<gcslpara2_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslpara2[i][0]=isurf;
        gcslpara2[i][1]=jsurf;
    }

    for(i=0; i<gcslpara3_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslpara3[i][0]=isurf;
        gcslpara3[i][1]=jsurf;
    }

    for(i=0; i<gcslpara4_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslpara4[i][0]=isurf;
        gcslpara4[i][1]=jsurf;
    }

    //  Paraslice Surfaces
    for(i=0; i<gcslparaco1_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslparaco1[i][0]=isurf;
        gcslparaco1[i][1]=jsurf;
    }

    for(i=0; i<gcslparaco2_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslparaco2[i][0]=isurf;
        gcslparaco2[i][1]=jsurf;
    }

    for(i=0; i<gcslparaco3_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslparaco3[i][0]=isurf;
        gcslparaco3[i][1]=jsurf;
    }

    for(i=0; i<gcslparaco4_count; ++i)
    {
        grid.read((char*)&iin, sizeof (int));
        isurf=iin;

        grid.read((char*)&iin, sizeof (int));
        jsurf=iin;

        gcslparaco4[i][0]=isurf;
        gcslparaco4[i][1]=jsurf;
    }

    // bed
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        grid.read((char*)&ddn, sizeof (double));
        bed[(i-imin)*jmax + (j-jmin)]=ddn;
    }

    if(solidread>0)
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        grid.read((char*)&ddn, sizeof (double));
        solidbed[(i-imin)*jmax + (j-jmin)]=ddn;
    }

    if(toporead>0)
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        grid.read((char*)&ddn, sizeof (double));
        topobed[(i-imin)*jmax + (j-jmin)]=ddn;
    }

    if(P150>0)
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
    {
        grid.read((char*)&ddn, sizeof (double));
        data[(i-imin)*jmax + (j-jmin)]=ddn;
    }

    #if USE_AMREX
    if(nlevs>1)
    {
        amrex_refined_grid_coords.resize(nlevs-1);
        for (int lev = 0; lev < nlevs-1; ++lev)
        {
            grid.read((char*)&iin, sizeof (int));
            int number_of_boxes = iin;
            amrex_refined_grid_coords[lev].resize(number_of_boxes);
            for (int box = 0; box < number_of_boxes; ++box)
            {
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].first[0] = ddn;
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].first[1] = ddn;
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].first[2] = ddn;
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].second[0] = ddn;
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].second[1] = ddn;
                grid.read((char*)&ddn, sizeof (double));
                amrex_refined_grid_coords[lev][box].second[2] = ddn;
            }
        }
        if(nlevs>max_nlevs) nlevs=max_nlevs;
    }
    #endif

    grid.close();
}
