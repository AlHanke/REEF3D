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

#include "printMethodCompact.h"

#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

#include "print_averaging.h"
#include "turbulence.h"
#include "heat.h"
#include "multiphase.h"
#include "vorticity.h"
#include "data.h"
#include "concentration.h"
#include "sediment.h"

#include <sys/stat.h>
#include <chrono>

printMethodCompact::printMethodCompact(lexer *p) : printMethod(p)
{
    if(p->mpirank==0)
        mkdir("./REEF3D_CFD_VTRC",0777);
}

printMethodCompact::~printMethodCompact()
{
    delete [] XN;
    delete [] YN;
    delete [] ZN;

    delete [] globalSendCounts;
    delete [] displs;

    delete [] uvel;
    delete [] vvel;
    delete [] wvel;
    delete [] press;
    delete [] eddyv;
    delete [] phi;
    delete [] topo;
    delete [] flag;
    delete [] flag4;
    delete [] flag5;

    delete [] pressGlobal;
    delete [] uvelGlobal;
    delete [] vvelGlobal;
    delete [] wvelGlobal;
    delete [] topoGlobal;
    delete [] phiGlobal;
    delete [] eddyvGlobal;
    delete [] flagGlobal;
    delete [] flag5Global;

}

void printMethodCompact::setup(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    int *recvcounts=nullptr;
    int *gneibours=nullptr;
    int *gextent=nullptr;

    if(p->mpirank==0)
    {
        recvcounts = new int[p->mpi_size];
        displs = new int[p->mpi_size];

        globalSendCounts = new int[p->mpi_size];

        gneibours = new int[p->mpi_size*6];
        gextent = new int[p->mpi_size*6];

        XN = new double[p->gknox+2];
        YN = new double[p->gknoy+2];
        ZN = new double[p->gknoz+2];

        cellNum = (p->gknox)*(p->gknoy)*(p->gknoz);
        pointNum = (p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);

        flag = new int*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        flag4 = new int*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        flag5 = new int*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];

        // ---------------------------------------------------------
        // Allocate memory for data to be printed
        // ---------------------------------------------------------
        
        uvel = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        vvel = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        wvel = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        press = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        eddyv = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        phi = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        topo = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P23==1)
            test = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P24==1 && p->F300==0)
            rho = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P25==1)
            solid = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P28==1)
            fb = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        // if(p->P29==1)
        //     walld = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P71==1)
            visc = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];
        if(p->P72==1)
            VOF = new double*[(p->gknox+2)*(p->gknoy+2)*(p->gknoz+2)];

        // ---------------------------------------------------------
        // Pre-calulate offsets
        // ---------------------------------------------------------

        n=0;
        vtkOffsets[n]=0;
        ++n;

        // calcVTKOffsets(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed,pointNum,cellNum);
        // Nothing besdides to following is implemented for printing

        //velocities
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*3*(pointNum);
        ++n;
        //pressure
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        //eddyv
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        //phi
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        // rho
        if(p->P24==1 && p->F300==0)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // viscosity
        if(p->P71==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // VOF
        if(p->P72==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // topo
        if(p->P27==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // test
        if(p->P23==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        //elevation
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        // solid
        if(p->P25==1)
        { 
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // floating
        if(p->P28==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        // walldist
        // if(p->P29==1)
        // {   
        //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        //     ++n;
        // }

        //x
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(p->gknox+1);
        ++n;
        //y
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(p->gknoy+1); 
        ++n;
        //z
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(p->gknoz+1);
    }

    int recvcount = p->knox+1;
    pgc->gather_int(&recvcount,1,recvcounts,1);
    int disp = p->origin_i+1;
    pgc->gather_int(&disp,1,displs,1);
    pgc->gatherv_double(p->XN+marge,p->knox+1, XN,recvcounts,displs);

    recvcount = p->knoy+1;
    pgc->gather_int(&recvcount,1,recvcounts,1);
    disp = p->origin_j+1;
    pgc->gather_int(&disp,1,displs,1);
    pgc->gatherv_double(p->YN+marge,p->knoy+1, YN,recvcounts,displs);

    recvcount = p->knoz+1;
    pgc->gather_int(&recvcount,1,recvcounts,1);
    disp = p->origin_k+1;
    pgc->gather_int(&disp,1,displs,1);
    pgc->gatherv_double(p->ZN+marge,p->knoz+1, ZN,recvcounts,displs);

    if(p->mpirank==0)
    {
        i=j=k=-1;
        XN[0]=p->XN[IP];
        YN[0]=p->YN[JP];
        ZN[0]=p->ZN[KP]; 
    }

    int neibours[6];
    neibours[0]=p->nb1;
    neibours[1]=p->nb2;
    neibours[2]=p->nb3;
    neibours[3]=p->nb4;
    neibours[4]=p->nb5;
    neibours[5]=p->nb6;
    pgc->gather_int(neibours,6,gneibours,6);

    int extent[6];
    extent[0]=p->origin_i;
    extent[1]=p->origin_i+p->knox;
    extent[2]=p->origin_j;
    extent[3]=p->origin_j+p->knoy;
    extent[4]=p->origin_k;
    extent[5]=p->origin_k+p->knoz;
    pgc->gather_int(extent,6,gextent,6);

    if(p->mpirank==0)
        localSendCount=0;
    else
        localSendCount=(p->knox+2*p->margin)*(p->knoy+2*p->margin)*(p->knoz+2*p->margin);
    pgc->gather_int(&localSendCount,1,globalSendCounts,1);

    if(p->mpirank==0)
    {
        int counter = globalSendCounts[0];
        displs[0]=0;
        for(int i=1;i<p->mpi_size;++i)
        {
            displs[i] = displs[i-1] + globalSendCounts[i-1];
            counter += globalSendCounts[i];
        }
        uvelGlobal = new double[counter];
        vvelGlobal = new double[counter];
        wvelGlobal = new double[counter];
        pressGlobal = new double[counter];
        eddyvGlobal = new double[counter];
        phiGlobal = new double[counter];
        topoGlobal = new double[counter];

        if(p->P23==1)
            testGlobal = new double[counter];
        if(p->P24==1 && p->F300==0)
            rhoGlobal = new double[counter];
        if(p->P25==1)
            solidGlobal = new double[counter];
        if(p->P28==1)
            fbGlobal = new double[counter];
        // if(p->P29==1)
        //     walldGlobal = new double[counter];
        if(p->P71==1)
            viscGlobal = new double[counter];
        if(p->P72==1)
            VOFGlobal = new double[counter];
        
        flagGlobal = new int[counter];
        flag4Global = new int[counter];
        flag5Global = new int[counter];
    }

    if(p->mpirank==0)
    {
        int indexL,indexLG;
        int kbegin,kend;
        int jbegin,jend;
        int ibegin,iend;
        for(int n=0;n<p->mpi_size;n++)
        {
            kbegin=-1;
            if(gneibours[4+6*n]>-2)
                kbegin=0;
            kend=gextent[5+6*n]-gextent[4+6*n]+1;
            if(gneibours[5+6*n]>-2)
                kend=gextent[5+6*n]-gextent[4+6*n];

            jbegin=-1;
            if(gneibours[2+6*n]>-2)
                jbegin=0;
            jend=gextent[3+6*n]-gextent[2+6*n]+1;
            if(gneibours[1+6*n]>-2)
                jend=gextent[3+6*n]-gextent[2+6*n];
            
            ibegin=-1;
            if(gneibours[0+6*n]>-2)
                ibegin=0;
            iend=gextent[1+6*n]-gextent[0+6*n]+1;
            if(gneibours[3+6*n]>-2)
                iend=gextent[1+6*n]-gextent[0+6*n];

            if(n!=0)
                for(int k=kbegin;k<kend;++k)
                {
                    for(int j=jbegin;j<jend;++j)
                    {
                        for(int i=ibegin;i<iend;++i)
                        {
                            indexL = (i+p->margin)*(gextent[3+6*n]-gextent[2+6*n]+2*p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + (j+p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + k+p->margin;
                            indexL += displs[n];
                            indexLG = (k+1+gextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+gextent[2+6*n])*(p->gknox+2)+(i+1+gextent[0+6*n]);
                            
                            uvel[indexLG] = &uvelGlobal[indexL];
                            vvel[indexLG] = &vvelGlobal[indexL];
                            wvel[indexLG] = &wvelGlobal[indexL];

                            press[indexLG] = &pressGlobal[indexL];
                            eddyv[indexLG] = &eddyvGlobal[indexL];
                            phi[indexLG] = &phiGlobal[indexL];
                            topo[indexLG] = &topoGlobal[indexL];
                                
                            if(p->P23==1)
                                test[indexLG] = &testGlobal[indexL];
                            if(p->P24==1 && p->F300==0)
                                rho[indexLG] = &rhoGlobal[indexL];
                            if(p->P25==1)
                                solid[indexLG] = &solidGlobal[indexL];
                            if(p->P28==1)
                                fb[indexLG] = &fbGlobal[indexL];
                            // if(p->P29==1)
                            //     walld[indexLG] = &walldGlobal[indexL];
                            if(p->P71==1)
                                visc[indexLG] = &viscGlobal[indexL];
                            if(p->P72==1)
                                VOF[indexLG] = &VOFGlobal[indexL];

                            flag[indexLG] = &flagGlobal[indexL];
                            flag4[indexLG] = &flag4Global[indexL];
                            flag5[indexLG] = &flag5Global[indexL];
                        }
                    }
                }
            else
                for(int k=kbegin;k<kend;++k)
                {
                    for(int j=jbegin;j<jend;++j)
                    {
                        for(int i=ibegin;i<iend;++i)
                        {
                            indexL = (i+p->margin)*(gextent[3+6*n]-gextent[2+6*n]+2*p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + (j+p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + k+p->margin;
                            indexLG = (k+1+gextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+gextent[2+6*n])*(p->gknox+2)+(i+1+gextent[0+6*n]);
                            
                            uvel[indexLG] = &a->u.V[indexL];
                            vvel[indexLG] = &a->v.V[indexL];
                            wvel[indexLG] = &a->w.V[indexL];
                            
                            press[indexLG] = &a->press.V[indexL];
                            eddyv[indexLG] = &a->eddyv.V[indexL];
                            phi[indexLG] = &a->phi.V[indexL];
                            topo[indexLG]=&a->topo.V[indexL];

                            if(p->P23==1)
                                test[indexLG]=&a->test.V[indexL];
                            if(p->P24==1 && p->F300==0)
                                rho[indexLG] = &a->ro.V[indexL];
                            if(p->P25==1)
                                solid[indexLG]=&a->solid.V[indexL];
                            if(p->P28==1)
                                fb[indexLG]=&a->fb.V[indexL];
                            // if(p->P29==1)
                            //     walld[indexLG]=&a->walld.V[indexL];
                            if(p->P71==1)
                                visc[indexLG] = &a->visc.V[indexL];
                            if(p->P72==1)
                                VOF[indexLG] = &a->vof.V[indexL];
                            
                            flag[indexLG]=&p->flag[indexL];
                            flag4[indexLG]=&p->flag4[indexL];
                            flag5[indexLG]=&p->flag5[indexL];
                        }
                    }
                }
        }
    }

    if(p->mpirank==0)
    {
        delete [] gneibours;
        delete [] gextent;
        delete [] recvcounts;
    }
}

int printMethodCompact::print(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    std::chrono::system_clock::time_point start,end;
    if(p->mpirank==0)
        start = std::chrono::system_clock::now();
    pgc->gatherv_double(a->u.V,localSendCount,uvelGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->v.V,localSendCount,vvelGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->w.V,localSendCount,wvelGlobal,globalSendCounts,displs);

    pgc->gatherv_double(a->press.V,localSendCount,pressGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->eddyv.V,localSendCount,eddyvGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->phi.V,localSendCount,phiGlobal,globalSendCounts,displs);

    if(p->P23==1)
        pgc->gatherv_double(a->test.V,localSendCount,testGlobal,globalSendCounts,displs);
    if(p->P24==1 && p->F300==0)
        pgc->gatherv_double(a->ro.V,localSendCount,rhoGlobal,globalSendCounts,displs);
    if(p->P25==1)
        pgc->gatherv_double(a->solid.V,localSendCount,solidGlobal,globalSendCounts,displs);
    // P27, also needed for others
    pgc->gatherv_double(a->topo.V,localSendCount,topoGlobal,globalSendCounts,displs);
    if(p->P28==1)
        pgc->gatherv_double(a->fb.V,localSendCount,fbGlobal,globalSendCounts,displs);
    // if(p->P29==1)
    //     pgc->gatherv_double(a->walld.V,localSendCount,walldGlobal,globalSendCounts,displs);
    if(p->P71==1)
        pgc->gatherv_double(a->visc.V,localSendCount,viscGlobal,globalSendCounts,displs);
    if(p->P72==1)
        pgc->gatherv_double(a->vof.V,localSendCount,VOFGlobal,globalSendCounts,displs);
    
    pgc->gatherv_int(p->flag,localSendCount,flagGlobal,globalSendCounts,displs);
    pgc->gatherv_int(p->flag4,localSendCount,flag4Global,globalSendCounts,displs);
    pgc->gatherv_int(p->flag5,localSendCount,flag5Global,globalSendCounts,displs);

    if(p->mpirank==0)
    {
        end = std::chrono::system_clock::now();
        auto elapsed = end - start;
        std::cout << "Communication time: "<<elapsed.count() << '\n';
    }

    int returnValue = 0;
    if(p->mpirank==0)
    {
        int num=0;
        if(p->P15==1)
            num = p->printcount;
        if(p->P15==2)
            num = p->count;
        sprintf(name,"./REEF3D_CFD_VTRC/REEF3D-CFD-%08i.vtr",num);
        n=0;
        ofstream result;
        start = std::chrono::system_clock::now();
        result.open(name,ios::binary);
        if(result.is_open())
        {
            result<<"<?xml version=\"1.0\"?>\n"
            <<"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
            <<"<RectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
            if(p->P16==1)
            {
            result<<"<FieldData>\n";
            result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<std::setprecision(7)<<p->simtime<<"</DataArray>\n";
            result<<"</FieldData>\n";
            }
            result<<"<Piece Extent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\">\n";
            result<<"<PointData>\n";
            result<<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"phi\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            if(p->P24==1 && p->F300==0)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"rho\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            if(p->P71==1)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"viscosity\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            if(p->P72==1)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"VOF\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            if(p->P27==1)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"topo\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            if(p->P23==1)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            result<<"\t<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            if(p->P25==1)
            { 
                result<<"\t<DataArray type=\"Float32\" Name=\"solid\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            if(p->P28==1)
            {
                result<<"\t<DataArray type=\"Float32\" Name=\"floating\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
                ++n;
            }
            // if(p->P29==1)
            // {   
            //     result<<"\t<DataArray type=\"Float32\" Name=\"walldist\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            //     ++n;
            // }
            result<<"</PointData>\n";
            result<<"<Coordinates>\n";
            result<<"\t<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\"/>\n";
            ++n;
            result<<"</Coordinates>\n"
            <<"</Piece>\n"
            <<"</RectilinearGrid>\n"
            <<"<AppendedData encoding=\"raw\">\n_";
            //  Velocities
            iin=3*4*(pointNum);
            result.write((char*)&iin, sizeof (int));
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol1(uvel,flag,flag5));//u
                        result.write((char*)&ffn, sizeof (float));

                        ffn=float(p->ipol2(vvel,flag,flag5));//v
                        result.write((char*)&ffn, sizeof (float));

                        ffn=float(p->ipol3(wvel,flag,flag5));//w
                        result.write((char*)&ffn, sizeof (float));
                    }
            //  Pressure
            iin=4*(pointNum);
            result.write((char*)&iin, sizeof (int));
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4press(press));
                        result.write((char*)&ffn, sizeof (float));
                    }
            //  EddyV
            iin=4*(pointNum);
            result.write((char*)&iin, sizeof (int));
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(eddyv));
                        result.write((char*)&ffn, sizeof (float));
                    }
            //  Phi
            iin=4*(pointNum);
            result.write((char*)&iin, sizeof (int));
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4phi(topo,phi));
                        result.write((char*)&ffn, sizeof (float));
                    }
            // rho
            if(p->P24==1 && p->F300==0)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4_a(rho));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // viscosity
            if(p->P71==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4(visc,flag4));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // VOF
            if(p->P72==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4(VOF,flag4));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // topo
            if(p->P27==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4_a(topo));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // test
            if(p->P23==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4_a(test));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            //  Elevation
            iin=4*(pointNum);
            result.write((char*)&iin, sizeof (int));
            for(k=0; k<p->gknoz+1; ++k)
                for(j=0; j<p->gknoy+1; ++j)
                    for(i=0; i<p->gknox+1; ++i)
                    {
                        ffn=float(ZN[k]+(ZN[k+1]-ZN[k]));
                        result.write((char*)&ffn, sizeof (float));
                    }
            // solid
            if(p->P25==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4_a(solid));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // floating
            if(p->P28==1)
            {
                iin=4*(pointNum);
                result.write((char*)&iin, sizeof (int));
                for(k=-1; k<p->gknoz; ++k)
                    for(j=-1; j<p->gknoy; ++j)
                        for(i=-1; i<p->gknox; ++i)
                        {
                            ffn=float(p->ipol4_a(fb));
                            result.write((char*)&ffn, sizeof (float));
                        }
            }
            // walldist
            // if(p->P29==1)
            // {
            //     iin=4*(pointNum);
            //     result.write((char*)&iin, sizeof (int));
            //     for(k=-1; k<p->gknoz; ++k)
            //         for(j=-1; j<p->gknoy; ++j)
            //             for(i=-1; i<p->gknox; ++i)
            //             {
            //                 ffn=float(p->ipol4_a(walld));
            //                 result.write((char*)&ffn, sizeof (float));
            //             }
            // }

            // x
            iin=4*(p->gknox+1);
            result.write((char*)&iin, sizeof (int));
            for(i=1; i<p->gknox+2; ++i)
            {
                ffn=float(XN[i]);
                result.write((char*)&ffn, sizeof (float));
            }
            // y
            iin=4*(p->gknoy+1);
            result.write((char*)&iin, sizeof (int));
            for(j=1; j<p->gknoy+2; ++j)
            {
                ffn=float(YN[j]);
                result.write((char*)&ffn, sizeof (float));
            }
            // z
            iin=4*(p->gknoz+1);
            result.write((char*)&iin, sizeof (int));
            for(k=1; k<p->gknoz+2; ++k)
            {
                ffn=float(ZN[k]);
                result.write((char*)&ffn, sizeof (float));
            }

            result<<"\n</AppendedData>\n"
            <<"</VTKFile>"<<flush;

            result.close();
            end = std::chrono::system_clock::now();
            auto elapsed = end - start;
            std::cout << "File print time time: "<<elapsed.count() << std::endl;
        }
        else
            returnValue = 1;
    }
    pgc->bcast_int(&returnValue,1);
    return returnValue;
}