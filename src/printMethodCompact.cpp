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

void printMethodCompact::setup(lexer *p, fdm *a, ghostcell *pgc)
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

        // ---------------------------------------------------------
        // Pre-calulate offsets
        // ---------------------------------------------------------

        n=0;
        vtkOffsets[n]=0;
        ++n;

        //velocities
        vtkOffsets[n]=vtkOffsets[n-1]+4+3*4*(pointNum);
        ++n;
        //pressure
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(pointNum);
        ++n;
        //eddyv
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(pointNum);
        ++n;
        //phi
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(pointNum);
        ++n;
        //elevation
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(pointNum);
        ++n;

        //x
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(p->gknox+1);
        ++n;
        //y
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(p->gknoy+1); 
        ++n;
        //z
        vtkOffsets[n]=vtkOffsets[n-1]+4+4*(p->gknoz+1);
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
        pressGlobal = new double[counter];
        uvelGlobal = new double[counter];
        vvelGlobal = new double[counter];
        wvelGlobal = new double[counter];
        topoGlobal = new double[counter];
        phiGlobal = new double[counter];
        eddyvGlobal = new double[counter];
        flagGlobal = new int[counter];
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
                            press[indexLG]=&pressGlobal[indexL];
                            uvel[indexLG]=&uvelGlobal[indexL];
                            vvel[indexLG]=&vvelGlobal[indexL];
                            wvel[indexLG]=&wvelGlobal[indexL];
                            topo[indexLG]=&topoGlobal[indexL];
                            phi[indexLG]=&phiGlobal[indexL];
                            eddyv[indexLG]=&eddyvGlobal[indexL];
                            flag[indexLG]=&flagGlobal[indexL];
                            flag5[indexLG]=&flag5Global[indexL];
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
                            press[indexLG]=&a->press.V[indexL];
                            uvel[indexLG]=&a->u.V[indexL];
                            vvel[indexLG]=&a->v.V[indexL];
                            wvel[indexLG]=&a->w.V[indexL];
                            topo[indexLG]=&a->topo.V[indexL];
                            phi[indexLG]=&a->phi.V[indexL];
                            eddyv[indexLG]=&a->eddyv.V[indexL];
                            flag[indexLG]=&p->flag[indexL];
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
    pgc->gatherv_double(a->press.V,localSendCount,pressGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->u.V,localSendCount,uvelGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->v.V,localSendCount,vvelGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->w.V,localSendCount,wvelGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->topo.V,localSendCount,topoGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->phi.V,localSendCount,phiGlobal,globalSendCounts,displs);
    pgc->gatherv_double(a->eddyv.V,localSendCount,eddyvGlobal,globalSendCounts,displs);
    pgc->gatherv_int(p->flag,localSendCount,flagGlobal,globalSendCounts,displs);
    pgc->gatherv_int(p->flag5,localSendCount,flag5Global,globalSendCounts,displs);

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
            result<<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            ++n;
            result<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            ++n;
            result<<"</PointData>\n";
            //  result<<"<CellData>\n";
            // <<"\t<DataArray type=\"Float32\" Name=\"Debug\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />\n";
            // ++n;
            // result<<"</CellData>\n";
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
        }
        else
            returnValue = 1;
    }
    pgc->bcast_int(&returnValue,1);
    return returnValue;
}