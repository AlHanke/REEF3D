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

#include <cstring>
#include <sstream>

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

    delete [] globalSendCountsFlag;
    delete [] globalSendCountsField;
    delete [] displs;
    delete [] displsField;
    delete [] displsFlag;

    delete [] domainSizes;
    delete [] allFlagsGlobal;
    delete [] allFlags;
    delete [] allFieldsGlobal;
    delete [] allFields;


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
    int *gindices=nullptr;
    int *gextent=nullptr;

    numberOfFlags=3;

    if(p->mpirank==0)
    {
        recvcounts = new int[p->mpi_size];
        displs = new int[p->mpi_size];
        displsField = new int[p->mpi_size];
        displsFlag = new int[p->mpi_size];

        globalSendCountsField = new int[p->mpi_size];
        globalSendCountsFlag = new int[p->mpi_size];

        gindices = new int[p->mpi_size*6];
        gextent = new int[p->mpi_size*6];
        domainSizes = new int[p->mpi_size];

        XN = new double[p->gknox+2];
        YN = new double[p->gknoy+2];
        ZN = new double[p->gknoz+2];

        cellNum = (p->gknox)*(p->gknoy)*(p->gknoz);
        pointNum = (p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);
        domainSize = (p->gknox+2)*(p->gknoy+2)*(p->gknoz+2);

        kbegin = p->nb5==-2?-1:0, kend = p->nb6==-2?p->knoz+1:p->knoz;
        jbegin = p->nb3==-2?-1:0, jend = p->nb2==-2?p->knoy+1:p->knoy;
        ibegin = p->nb1==-2?-1:0, iend = p->nb4==-2?p->knox+1:p->knox;

        flag = new int*[domainSize];
        flag4 = new int*[domainSize];
        flag5 = new int*[domainSize];

        // ---------------------------------------------------------
        // Allocate memory for data to be printed
        // ---------------------------------------------------------
        
        uvel = new double*[domainSize];
        vvel = new double*[domainSize];
        wvel = new double*[domainSize];
        press = new double*[domainSize];
        eddyv = new double*[domainSize];
        phi = new double*[domainSize];
        topo = new double*[domainSize];
        if(p->P23==1)
            test = new double*[domainSize];
        if(p->P24==1 && p->F300==0)
            rho = new double*[domainSize];
        if(p->P25==1)
            solid = new double*[domainSize];
        if(p->P28==1)
            fb = new double*[domainSize];
        // if(p->P29==1)
        //     walld = new double*[domainSize];
        if(p->P71==1)
            visc = new double*[domainSize];
        if(p->P72==1)
            VOF = new double*[domainSize];

        // ---------------------------------------------------------
        // Pre-calulate offsets
        // ---------------------------------------------------------

        n=0;
        vtkOffsets[n]=0;
        ++n;
        numberOfFields = 0;

        // calcVTKOffsets(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed,pointNum,cellNum);
        // Nothing besdides to following is implemented for printing

        //velocities
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*3*(pointNum);
        ++n;
        numberOfFields += 3;
        //pressure
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        numberOfFields += 1;
        //eddyv
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        numberOfFields += 1;
        //phi
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        numberOfFields += 1;
        // rho
        if(p->P24==1 && p->F300==0)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        // viscosity
        if(p->P71==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        // VOF
        if(p->P72==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        // topo
        if(p->P27==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
        }
        numberOfFields += 1;
        // test
        if(p->P23==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        //elevation
        vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++n;
        // solid
        if(p->P25==1)
        { 
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        // floating
        if(p->P28==1)
        {
            vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
            ++n;
            numberOfFields += 1;
        }
        // walldist
        // if(p->P29==1)
        // {   
        //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(float)*(pointNum);
        //     ++n;
        //     numberOfFields += 1;
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
    else
    {
        cellNum = p->knox*p->knoy*p->knoz;
        pointNum = (p->knox+1)*(p->knoy+1)*(p->knoz+1);
        kbegin = p->nb5==-2?-1:0, kend = p->nb6==-2?p->knoz+1:p->knoz;
        jbegin = p->nb3==-2?-1:0, jend = p->nb2==-2?p->knoy+1:p->knoy;
        ibegin = p->nb1==-2?-1:0, iend = p->nb4==-2?p->knox+1:p->knox;
        domainSize = (iend-ibegin)*(jend-jbegin)*(kend-kbegin);
    }

    pgc->bcast_int(&numberOfFields,1);

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

    pgc->gather_int(&domainSize,1,domainSizes,1);

    if(p->mpirank==0)
    {
        i=j=k=-1;
        XN[0]=p->XN[IP];
        YN[0]=p->YN[JP];
        ZN[0]=p->ZN[KP]; 
    }

    if(p->mpirank==0)
    {
        domainSize -= (iend-ibegin)*(jend-jbegin)*(kend-kbegin);
        allFieldsGlobal = new double [numberOfFields*domainSize];
        allFields = new double* [numberOfFields*domainSize];
        allFlagsGlobal = new int [numberOfFlags*domainSize];
    }
    else
    {
        allFields = new double* [numberOfFields*domainSize];
        allFieldsGlobal = new double [numberOfFields*domainSize];
        allFlags = new int*[numberOfFlags*domainSize];
        allFlagsGlobal = new int [numberOfFlags*domainSize];
        
        int m;
        n=0;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->u(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->v(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->w(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->press(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->eddyv(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->phi(i,j,k);
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFields[m+n*domainSize] = &a->topo(i,j,k);
            ++m;
        }
        ++n;

        n=0;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFlags[m+n*domainSize] = &p->flag[IJK];
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFlags[m+n*domainSize] = &p->flag4[IJK];
            ++m;
        }
        ++n;
        m=0;
        for(k=kbegin; k<kend; ++k)
        for(j=jbegin; j<jend; ++j)
        for(i=ibegin; i<iend; ++i)
        {
            allFlags[m+n*domainSize] = &p->flag5[IJK];
            ++m;
        }
        ++n;
    }

    int indices[6];
    indices[0]=kbegin;
    indices[1]=kend;
    indices[2]=jbegin;
    indices[3]=jend;
    indices[4]=ibegin;
    indices[5]=iend;
    pgc->gather_int(indices,6,gindices,6);

    int extent[6];
    extent[0]=p->origin_i;
    extent[1]=p->origin_i+p->knox;
    extent[2]=p->origin_j;
    extent[3]=p->origin_j+p->knoy;
    extent[4]=p->origin_k;
    extent[5]=p->origin_k+p->knoz;
    pgc->gather_int(extent,6,gextent,6);

    if(p->mpirank==0)
    {
        localSendCountField=0;
        localSendCountFlag=0;
    }
    else
    {
        localSendCountField = numberOfFields * domainSize;
        localSendCountFlag = numberOfFlags * domainSize;
    }
    pgc->gather_int(&localSendCountField,1,globalSendCountsField,1);
    pgc->gather_int(&localSendCountFlag,1,globalSendCountsFlag,1);

    if(p->mpirank==0)
    {
        displsField[0]=0;
        displsFlag[0]=0;
        for(int i=1;i<p->mpi_size;++i)
        {
            displsField[i] = displsField[i-1] + globalSendCountsField[i-1];
            displsFlag[i] = displsFlag[i-1] + globalSendCountsFlag[i-1];
        }
    }

    if(p->mpirank==0)
    {
        int indexL,indexLG;
        for(int n=0;n<p->mpi_size;n++)
        {
            for(int k=gindices[0+6*n];k<gindices[1+6*n];++k)
            {
                for(int j=gindices[2+6*n];j<gindices[3+6*n];++j)
                {
                    for(int i=gindices[4+6*n];i<gindices[5+6*n];++i)
                    {
                        if(n!=0)
                        {
                            int m=0;
                            indexL = (i+p->margin)*(gextent[3+6*n]-gextent[2+6*n]+2*p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + (j+p->margin)*(gextent[5+6*n]-gextent[4+6*n]+2*p->margin) + k+p->margin;
                            indexL += displs[n];
                            indexLG = (k+1+gextent[4+6*n])*(p->gknox+2)*(p->gknoy+2)+(j+1+gextent[2+6*n])*(p->gknox+2)+(i+1+gextent[0+6*n]);
                            
                            uvel[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            vvel[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            wvel[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;

                            press[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            eddyv[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            phi[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            topo[indexLG] = &allFieldsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                                
                            // if(p->P23==1)
                            //     test[indexLG] = &testGlobal[indexL];
                            // if(p->P24==1 && p->F300==0)
                            //     rho[indexLG] = &rhoGlobal[indexL];
                            // if(p->P25==1)
                            //     solid[indexLG] = &solidGlobal[indexL];
                            // if(p->P28==1)
                            //     fb[indexLG] = &fbGlobal[indexL];
                            // // if(p->P29==1)
                            // //     walld[indexLG] = &walldGlobal[indexL];
                            // if(p->P71==1)
                            //     visc[indexLG] = &viscGlobal[indexL];
                            // if(p->P72==1)
                            //     VOF[indexLG] = &VOFGlobal[indexL];

                            m=0;
                            flag[indexLG] = &allFlagsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            flag4[indexLG] = &allFlagsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                            flag5[indexLG] = &allFlagsGlobal[indexL+domainSizes[n]*m];
                            ++m;
                        }
                        else
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

                            // if(p->P23==1)
                            //     test[indexLG]=&a->test.V[indexL];
                            // if(p->P24==1 && p->F300==0)
                            //     rho[indexLG] = &a->ro.V[indexL];
                            // if(p->P25==1)
                            //     solid[indexLG]=&a->solid.V[indexL];
                            // if(p->P28==1)
                            //     fb[indexLG]=&a->fb.V[indexL];
                            // // if(p->P29==1)
                            // //     walld[indexLG]=&a->walld.V[indexL];
                            // if(p->P71==1)
                            //     visc[indexLG] = &a->visc.V[indexL];
                            // if(p->P72==1)
                            //     VOF[indexLG] = &a->vof.V[indexL];
                            
                            flag[indexLG]=&p->flag[indexL];
                            flag4[indexLG]=&p->flag4[indexL];
                            flag5[indexLG]=&p->flag5[indexL];
                        }
                    }
                }
            }
        }
    }

    if(p->mpirank==0)
    {
        delete [] gindices;
        delete [] gextent;
        delete [] recvcounts;
    }
}

void printMethodCompact::fillContainers(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    for (int n=0;n<numberOfFields*domainSize;n++)
        allFieldsGlobal[n] = *allFields[n];
    for (int n=0;n<numberOfFlags*domainSize;n++)
        allFlagsGlobal[n] = *allFlags[n];
}

int printMethodCompact::print(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    int returnValue = 0;
    std::chrono::system_clock::time_point start,end;

    // pgc->gatherv_double(a->u.V,localSendCount,uvelGlobal,globalSendCounts,displs);
    // pgc->gatherv_double(a->v.V,localSendCount,vvelGlobal,globalSendCounts,displs);
    // pgc->gatherv_double(a->w.V,localSendCount,wvelGlobal,globalSendCounts,displs);

    // pgc->gatherv_double(a->press.V,localSendCount,pressGlobal,globalSendCounts,displs);
    // pgc->gatherv_double(a->eddyv.V,localSendCount,eddyvGlobal,globalSendCounts,displs);
    // pgc->gatherv_double(a->phi.V,localSendCount,phiGlobal,globalSendCounts,displs);

    // if(p->P23==1)
    //     pgc->gatherv_double(a->test.V,localSendCount,testGlobal,globalSendCounts,displs);
    // if(p->P24==1 && p->F300==0)
    //     pgc->gatherv_double(a->ro.V,localSendCount,rhoGlobal,globalSendCounts,displs);
    // if(p->P25==1)
    //     pgc->gatherv_double(a->solid.V,localSendCount,solidGlobal,globalSendCounts,displs);
    // // P27, also needed for others
    // pgc->gatherv_double(a->topo.V,localSendCount,topoGlobal,globalSendCounts,displs);
    // if(p->P28==1)
    //     pgc->gatherv_double(a->fb.V,localSendCount,fbGlobal,globalSendCounts,displs);
    // // if(p->P29==1)
    // //     pgc->gatherv_double(a->walld.V,localSendCount,walldGlobal,globalSendCounts,displs);
    // if(p->P71==1)
    //     pgc->gatherv_double(a->visc.V,localSendCount,viscGlobal,globalSendCounts,displs);
    // if(p->P72==1)
    //     pgc->gatherv_double(a->vof.V,localSendCount,VOFGlobal,globalSendCounts,displs);

    
    
    // pgc->gatherv_int(p->flag,localSendCount,flagGlobal,globalSendCounts,displs);
    // pgc->gatherv_int(p->flag4,localSendCount,flag4Global,globalSendCounts,displs);
    // pgc->gatherv_int(p->flag5,localSendCount,flag5Global,globalSendCounts,displs);

    
    if(p->mpirank==0)
    {
        start = std::chrono::system_clock::now();
        std::stringstream result;
        n=0;
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

        m=result.str().length();
        buffer.resize(m+vtkOffsets[n]+27);
        std::memcpy(&buffer[0],result.str().data(),m);
        end = std::chrono::system_clock::now();
        auto elapsed = end - start;
        std::cout << "Header time: "<<elapsed.count() << '\n';
    }
    else
    {
        start = std::chrono::system_clock::now();
        fillContainers(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed);
        end = std::chrono::system_clock::now();
        auto elapsed = end - start;
        std::cout << "Fill container time of "<<p->mpirank<<": "<<elapsed.count() << '\n';
    }

    if(p->mpirank==0)
        start = std::chrono::system_clock::now();
    pgc->gatherv_double(allFieldsGlobal,localSendCountField,allFieldsGlobal,globalSendCountsField,displsField);
    pgc->gatherv_int(allFlagsGlobal,localSendCountFlag,allFlagsGlobal,globalSendCountsFlag,displsFlag);

    if(p->mpirank==0)
    {
        end = std::chrono::system_clock::now();
        auto elapsed = end - start;
        std::cout << "Communication time: "<<elapsed.count() << '\n';
    }

    
    if(p->mpirank==0)
    {
        start = std::chrono::system_clock::now();
        //  Velocities
        iin=3*4*(pointNum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=-1; k<p->gknoz; ++k)
            for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol1(uvel,flag,flag5));//u
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);

                    ffn=float(p->ipol2(vvel,flag,flag5));//v
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);

                    ffn=float(p->ipol3(wvel,flag,flag5));//w
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);
                }
        //  Pressure
        iin=4*(pointNum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=-1; k<p->gknoz; ++k)
            for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol4press(press));
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);
                }
        //  EddyV
        iin=4*(pointNum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=-1; k<p->gknoz; ++k)
            for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol4_a(eddyv));
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);
                }
        //  Phi
        iin=4*(pointNum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=-1; k<p->gknoz; ++k)
            for(j=-1; j<p->gknoy; ++j)
                for(i=-1; i<p->gknox; ++i)
                {
                    ffn=float(p->ipol4phi(topo,phi));
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);
                }
        // rho
        if(p->P24==1 && p->F300==0)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(rho));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // viscosity
        if(p->P71==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4(visc,flag4));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // VOF
        if(p->P72==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4(VOF,flag4));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // topo
        if(p->P27==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(topo));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // test
        if(p->P23==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(test));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        //  Elevation
        iin=4*(pointNum);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=0; k<p->gknoz+1; ++k)
            for(j=0; j<p->gknoy+1; ++j)
                for(i=0; i<p->gknox+1; ++i)
                {
                    ffn=float(ZN[k]+(ZN[k+1]-ZN[k]));
                    std::memcpy(&buffer[m],&ffn,sizeof(float));
                    m+=sizeof(float);
                }
        // solid
        if(p->P25==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(solid));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // floating
        if(p->P28==1)
        {
            iin=4*(pointNum);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            for(k=-1; k<p->gknoz; ++k)
                for(j=-1; j<p->gknoy; ++j)
                    for(i=-1; i<p->gknox; ++i)
                    {
                        ffn=float(p->ipol4_a(fb));
                        std::memcpy(&buffer[m],&ffn,sizeof(float));
                        m+=sizeof(float);
                    }
        }
        // walldist
        // if(p->P29==1)
        // {
        //     iin=4*(pointNum);
        //     std::memcpy(&buffer[m],&iin,sizeof(int));
        //     m+=sizeof(int);
        //     for(k=-1; k<p->gknoz; ++k)
        //         for(j=-1; j<p->gknoy; ++j)
        //             for(i=-1; i<p->gknox; ++i)
        //             {
        //                 ffn=float(p->ipol4_a(walld));
        //                 std::memcpy(&buffer[m],&ffn,sizeof(float));
        //                 m+=sizeof(float);
        //             }
        // }

        // x
        iin=4*(p->gknox+1);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(i=1; i<p->gknox+2; ++i)
        {
            ffn=float(XN[i]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }
        // y
        iin=4*(p->gknoy+1);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(j=1; j<p->gknoy+2; ++j)
        {
            ffn=float(YN[j]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }
        // z
        iin=4*(p->gknoz+1);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(k=1; k<p->gknoz+2; ++k)
        {
            ffn=float(ZN[k]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }
        std::memcpy(&buffer[m],&"\n</AppendedData>\n</VTKFile>",28);
        end = std::chrono::system_clock::now();
        auto elapsed = end - start;
        std::cout << "Content time: "<<elapsed.count() << std::endl;

        int num=0;
        if(p->P15==1)
            num = p->printcount;
        if(p->P15==2)
            num = p->count;
        sprintf(name,"./REEF3D_CFD_VTRC/REEF3D-CFD-%08i.vtr",num);
        start = std::chrono::system_clock::now();
        FILE* file = fopen(name, "w");
        fwrite(buffer.data(), buffer.size(), 1, file);
        fclose(file);
        end = std::chrono::system_clock::now();
        elapsed = end - start;
        std::cout << "File print time: "<<elapsed.count() << std::endl;
    }
    pgc->bcast_int(&returnValue,1);
    return returnValue;
}