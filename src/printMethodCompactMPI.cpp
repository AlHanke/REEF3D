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

#include "printMethodCompactMPI.h"

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

printMethodCompactMPI::printMethodCompactMPI(lexer *p) : printMethod(p), eta(p), node(p)
{
    if(p->mpirank==0)
        mkdir("./REEF3D_CFD_VTRCMPI",0777);
}

void printMethodCompactMPI::setup(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    int *gorigins=nullptr;
    int *gbeginEndPoint = nullptr;

    double *XN = nullptr;
    double *YN = nullptr;
    double *ZN = nullptr;

    int *recvcounts=nullptr;
    int *displs=nullptr;


    if(p->mpirank==0)
    {
	    mkdir("./REEF3D_CFD_VTRCMPI",0777);

        XN = new double[p->gknox+2];
        YN = new double[p->gknoy+2];
        ZN = new double[p->gknoz+2];

        gbeginEndPoint = new int[p->mpi_size*6];
        gorigins = new int[p->mpi_size*3];

        recvcounts = new int[p->mpi_size];
        displs = new int[p->mpi_size];

        cellNum=(p->gknox)*(p->gknoy)*(p->gknoz);
        pointNum=(p->gknox+1)*(p->gknoy+1)*(p->gknoz+1);

        // ---------------------------------------------------------
        // Pre-calulate offsets
        // ---------------------------------------------------------

        m=0;
        compactMPIPOffset[m]=0;
        ++m;

        // time
        if(p->P16==1)
        {
            compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(double);
            ++m;
        }

        //velocities
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+3*sizeof(float)*(pointNum);
        ++m;
        //pressure
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++m;
        //eddyv
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++m;
        //phi
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++m;
        //elevation
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(pointNum);
        ++m;

        //x
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknox+1);
        ++m;
        //y
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknoy+1); 
        ++m;
        //z
        compactMPIPOffset[m]=compactMPIPOffset[m-1]+sizeof(int)+sizeof(float)*(p->gknoz+1);
        ++m;
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


    kbeginPoint=-1;
    kendPoint=p->knoz;
    if(p->nb6>-2)
        --kendPoint;

    jbeginPoint=-1;
    jendPoint=p->knoy;
    if(p->nb2>-2)
        --jendPoint;
    
    ibeginPoint=-1;
    iendPoint=p->knox;
    if(p->nb4>-2)
        --iendPoint;

    int beginEndPoint[6];
    beginEndPoint[0]=ibeginPoint;
    beginEndPoint[1]=iendPoint;
    beginEndPoint[2]=jbeginPoint;
    beginEndPoint[3]=jendPoint;
    beginEndPoint[4]=kbeginPoint;
    beginEndPoint[5]=kendPoint;
    pgc->gather_int(beginEndPoint,6,gbeginEndPoint,6);
    int origins[3];
    origins[0]=p->origin_i;
    origins[1]=p->origin_j;
    origins[2]=p->origin_k;
    pgc->gather_int(origins,3,gorigins,3);

    if(p->mpirank==0)
    {
        // ---------------------------------------------------------
        // Data MPI offsets
        // ---------------------------------------------------------
        {
            m=0;

            //time
            if(p->P16==1)
            {
                offsetCMPI.push_back(compactMPIPOffset[m]);
                for(int n=0;n<p->mpi_size;++n)
                {
                    offsetCMPIitr.push_back(offsetCMPI.size()-1);
                    offsetCMPI.push_back(compactMPIPOffset[m]);
                }
                offsetCMPI.pop_back();
                m++;
            }

            //velocities
            offsetCMPIPoints(p,gorigins,gbeginEndPoint,3);

            //pressure
            offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

            //eddyv
            offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

            //phi
            offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);

            //elevation
            offsetCMPIPoints(p,gorigins,gbeginEndPoint,1);
        }

        // header
        {
            header<<"<?xml version=\"1.0\"?>\n"
            <<"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
            <<"<RectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
            m=0;
            if(p->P16==1)
            {
                header<<"\t<FieldData>\n";
                header<<"\t\t<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
                ++m;
                header<<"\t</FieldData>\n";
            }
            header<<"\t<Piece Extent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\">\n";
            header<<"\t\t<PointData>\n";
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
            ++m;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
            ++m;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
            ++m;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"phi\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
            ++m;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\" />\n";
            ++m;
            header<<"\t\t</PointData>\n";
            endIndex=m;
            header<<"\t\t<Coordinates>\n";
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
            m++;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
            m++;
            header<<"\t\t\t<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<compactMPIPOffset[m]<<"\"/>\n";
            m++;
            header<<"\t\t</Coordinates>\n"
            <<"\t</Piece>\n"
            <<"</RectilinearGrid>\n"
            <<"<AppendedData encoding=\"raw\">\n_";
            headerSize=header.str().size();
            for(auto &elem : offsetCMPI)
                elem+=headerSize;
        }

        // footer
        {
            // x
            iin=4*(p->gknox+1);
            footer.write((char*)&iin, sizeof (int));
            for(i=1; i<p->gknox+2; ++i)
            {
                ffn=float(XN[i]);
                footer.write((char*)&ffn, sizeof (float));
            }
            // y
            iin=4*(p->gknoy+1);
            footer.write((char*)&iin, sizeof (int));
            for(j=1; j<p->gknoy+2; ++j)
            {
                ffn=float(YN[j]);
                footer.write((char*)&ffn, sizeof (float));
            }
            // z
            iin=4*(p->gknoz+1);
            footer.write((char*)&iin, sizeof (int));
            for(k=1; k<p->gknoz+2; ++k)
            {
                ffn=float(ZN[k]);
                footer.write((char*)&ffn, sizeof (float));
            }

            footer<<"\n</AppendedData>\n"
            <<"</VTKFile>";
        }
    }
    
    int size=0;
    if(p->mpirank==0)
        size = offsetCMPI.size();
    pgc->Bcast(&size,1,MPI_INT);
    if(p->mpirank!=0)
        offsetCMPI.resize(size);
    pgc->Bcast(offsetCMPI.data(),size,MPI_OFFSET);

    if(p->mpirank==0)
        size = offsetCMPIitr.size();
    pgc->Bcast(&size,1,MPI_INT);
    if(p->mpirank!=0)
        offsetCMPIitr.resize(size);
    pgc->Bcast(offsetCMPIitr.data(),size,MPI_INT);

    if(p->mpirank==0)
    {
        delete[] XN;
        delete[] YN;
        delete[] ZN;

        delete[] gorigins;
        delete[] gbeginEndPoint;
        delete[] recvcounts;
        delete[] displs;
    }
}

void printMethodCompactMPI::offsetCMPIPoints(lexer *p, int *gorigins, int *gbeginEndPoint, int numberOfTuples)
{
    offsetCMPI.push_back(compactMPIPOffset[m]);
    for(int n=0;n<p->mpi_size;++n)
    {
        offsetCMPIitr.push_back(offsetCMPI.size()-1);
        if(n>0)
        {
            offsetCMPI.pop_back();
            offsetCMPI.push_back(compactMPIPOffset[m]+(gorigins[0+3*n]+gorigins[1+3*n]*(p->gknox+1)+gorigins[2+3*n]*(p->gknox+1)*(p->gknoy+1))*numberOfTuples*sizeof(float)+sizeof(int));
        }
        for(int k=0;k<(gbeginEndPoint[5+6*n]-gbeginEndPoint[4+6*n]);k++)
        {
            if(k>0)
            {
                offsetCMPI.pop_back();
                offsetCMPI.push_back(compactMPIPOffset[m]+(gorigins[0+3*n]+gorigins[1+3*n]*(p->gknox+1)+(gorigins[2+3*n]+k)*(p->gknox+1)*(p->gknoy+1))*numberOfTuples*sizeof(float)+sizeof(int));
            }
            for(int j=0;j<(gbeginEndPoint[3+6*n]-gbeginEndPoint[2+6*n]);j++)
            {
                offsetCMPI.push_back(offsetCMPI.back()+numberOfTuples*sizeof(float)*(p->gknox+1));
                if(n==0&&k==0&&j==0)
                    offsetCMPI.back()+=sizeof(int);
            }
        }
    }
    offsetCMPI.pop_back();
    m++;
}

int printMethodCompactMPI::print(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    pgc->start4a(p,a->test,1);
    pgc->start1(p,a->u,110);
    pgc->start2(p,a->v,111);
    pgc->start3(p,a->w,112);


    pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
    pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
    pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol4(p,a->phi,14);
    pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
    //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

    a->u.ggcpol(p);
    a->v.ggcpol(p);
    a->w.ggcpol(p);
    a->press.ggcpol(p);
    a->eddyv.ggcpol(p);
    a->phi.ggcpol(p);
    a->conc.ggcpol(p);
    a->ro.ggcpol(p);
    a->visc.ggcpol(p);
    a->phi.ggcpol(p);
    a->fb.ggcpol(p);
    a->fbh4.ggcpol(p);
    //a->test.ggcpol(p);
    

    pgc->gcparacox(p,a->phi,50);
    pgc->gcparacox(p,a->phi,50);

    pgc->gcparacox(p,a->topo,150);
    pgc->gcparacox(p,a->topo,150);
    
    //pgc->start4a(p,a->topo,159);

    pgc->gcperiodicx(p,a->press,4);


    MPI_File file;
    if(p->mpirank==0)
    {
        int num=0;
        if(p->P15==1)
            num = p->printcount;
        if(p->P15==2)
            num = p->count;
        snprintf(name,200,"./REEF3D_CFD_VTRCMPI/REEF3D-CFD-%08i.vtr",num);
    }
    pgc->Bcast(name,200,MPI_CHAR);
    pgc->File_open_createWriteOnly(&file,name);
    pgc->File_set_size(file,0);

    // header
    if(p->mpirank==0)
    {
        pgc->File_write_at_char(file, 0, header.str().c_str(), header.str().size());
    }

    // data + footer
    {
        std::stringstream result;
        m=0;

        // Time
        if(p->P16==1)
        {
            if(p->mpirank==0)
            {
                iin=sizeof(double);
                result.write((char*)&iin, sizeof (int));
                std::stringstream time;
                time<<std::setprecision(7)<<p->simtime;
                double t = std::stod(time.str());
                result.write((char*)&t, sizeof(double));
                pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m]], result.str().c_str(), result.str().size());
                result.str(std::string());
                result.clear();
            }
            ++m;
        }

        //  Velocities
        {
            if(p->mpirank==0)
            {
                iin=3*sizeof(float)*(pointNum);
                result.write((char*)&iin, sizeof (int));
            }
            int n=0;
            for(k=kbeginPoint;k<kendPoint;++k)
                for(j=jbeginPoint;j<jendPoint;++j)
                {
                    for(i=ibeginPoint;i<iendPoint;++i)
                    {
                        ffn=float(p->ipol1(a->u));
                        result.write((char*)&ffn, sizeof (float));

                        ffn=float(p->ipol2(a->v));
                        result.write((char*)&ffn, sizeof (float));

                        ffn=float(p->ipol3(a->w));
                        result.write((char*)&ffn, sizeof (float));
                    }
                    pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
                    result.str(std::string());
                    result.clear();
                    ++n;
                }
            ++m;
        }

        //  Pressure
        {
            if(p->mpirank==0)
            {
                iin=3*sizeof(float)*(pointNum);
                result.write((char*)&iin, sizeof (int));
            }
            int n=0;
            for(k=kbeginPoint;k<kendPoint;++k)
                for(j=jbeginPoint;j<jendPoint;++j)
                {
                    for(i=ibeginPoint;i<iendPoint;++i)
                    {
                        ffn=float(p->ipol4press(a->press)-p->pressgage);
                        result.write((char*)&ffn, sizeof (float));
                    }
                    pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
                    result.str(std::string());
                    result.clear();
                    ++n;
                }
            ++m;
        }

        //  EddyV
        {
            if(p->mpirank==0)
            {
                iin=sizeof(float)*(pointNum);
                result.write((char*)&iin, sizeof (int));
            }
            int n=0;
            for(k=kbeginPoint;k<kendPoint;++k)
                for(j=jbeginPoint;j<jendPoint;++j)
                {
                    for(i=ibeginPoint;i<iendPoint;++i)
                    {
                        ffn=float(p->ipol4_a(a->eddyv));
                        result.write((char*)&ffn, sizeof (float));
                    }
                    pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
                    result.str(std::string());
                    result.clear();
                    ++n;
                }
            ++m;
        }

        //  Phi
        {
            if(p->mpirank==0)
            {
                iin=sizeof(float)*(pointNum);
                result.write((char*)&iin, sizeof (int));
            }
            int n=0;
            node.nodefill4(p,a,pgc,a->phi,eta);
            for(k=kbeginPoint;k<kendPoint;++k)
                for(j=jbeginPoint;j<jendPoint;++j)
                {
                    for(i=ibeginPoint;i<iendPoint;++i)
                    {
                        if(p->P18==1)
                            ffn=float(p->ipol4phi(a,a->phi));
                        if(p->P18==2)
                            ffn = float(eta(i,j,k));
                        result.write((char*)&ffn, sizeof (float));
                    }
                    pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
                    result.str(std::string());
                    result.clear();
                    ++n;
                }
            ++m;
        }

        //  Elevation
        {
            if(p->mpirank==0)
            {
                iin=sizeof(float)*(pointNum);
                result.write((char*)&iin, sizeof (int));
            }
            int n=0;
            for(k=kbeginPoint;k<kendPoint;++k)
                for(j=jbeginPoint;j<jendPoint;++j)
                {
                    for(i=ibeginPoint;i<iendPoint;++i)
                    {
                        ffn=float(p->pos_z()+0.5*p->DZN[KP]);
                        result.write((char*)&ffn, sizeof (float));
                    }
                    pgc->File_write_at_char(file, offsetCMPI[offsetCMPIitr[m*p->mpi_size+p->mpirank]+n], result.str().c_str(), result.str().size());
                    result.str(std::string());
                    result.clear();
                    ++n;
                }
            ++m;
        }

        // footer
        if(p->mpirank==0)
        {
            pgc->File_write_at_char(file, headerSize + compactMPIPOffset[m], footer.str().c_str(), footer.str().size());
        }
    }
    pgc->File_close(&file);
    return 0;
}