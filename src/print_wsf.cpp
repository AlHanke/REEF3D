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

#include "print_wsf.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

print_wsf::print_wsf(lexer *p, fdm* a, ghostcell *pgc, int num)
{
    if(p->P51>0 && num==0)
    {
        gauge_num = p->P51;
        x = p->P51_x;
        y = p->P51_y;
    }
    else if(p->P351>0 && num==1)
    {
        gauge_num = p->P351;
        x = p->P351_x;
        y = p->P351_y;
    }
    else if(p->P352>0 && num==2)
    {
        gauge_num = p->P352;
        x = p->P352_x;
        y = p->P352_y;
    }
    else
    {
        if(p->mpirank==0)
        std::cerr<<"Error: Impropper height gauge defined!"<<std::endl;

        pgc->final(EXIT_FAILURE);
    }

    if(p->mpirank==0)
    {
        // Create Folder
        mkdir("./REEF3D_CFD_WSF",0777);

        // open file
        if(num==0)
        wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG.dat");

        else if(num==1)
        wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-1.dat");

        else if(num==2)
        wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-2.dat");
    }

    if(p->mpirank==0)
    {
        wsfout<<"number of gauges:  "<<gauge_num<<"\n\n";
        wsfout<<"x_coord     y_coord\n";
        for(n=0;n<gauge_num;++n)
        wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<"\n\n\n";

        wsfout<<"time";
        for(n=0;n<gauge_num;++n)
        wsfout<<"\t P"<<n+1;

        wsfout<<"\n\n"<<std::flush;
    }

    iloc.resize(gauge_num);
    jloc.resize(gauge_num);
    flag.resize(gauge_num);
    wsf.resize(gauge_num);

    ini_location(p);
}

print_wsf::~print_wsf()
{
    wsfout.close();
}

void print_wsf::height_gauge(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    std::fill(wsf.begin(), wsf.end(), -1.0e20);

    if(p->A10==6 && p->F80!=4)
    {
        for(n=0;n<gauge_num;++n)
        if(flag[n]>0)
        {
            i=iloc[n];
            j=jloc[n];
            KLOOP
            PCHECK
            {
                if(f(i,j,k)>=0.0 && f(i,j,k+1)<0.0)
                wsf[n]=MAX(wsf[n],-(f(i,j,k)*p->DZP[KP])/(f(i,j,k+1)-f(i,j,k)) + p->pos_z());
            }
        }
    }
    else if(p->A10==6 && p->F80==4)
    {
        for(n=0;n<gauge_num;++n)
        if(flag[n]>0)
        {
            i=iloc[n];
            j=jloc[n];
            KLOOP
            {
                if(f(i,j,k)>p->F94 && f(i,j,k+1)<p->F93)
                wsf[n]=MAX(wsf[n],p->pos_z()+0.5*p->DZN[KP]);

                else if(f(i,j,k)<=p->F94 && f(i,j,k)>=p->F93)
                wsf[n]=MAX(wsf[n],(p->pos_z()-0.5*p->DZN[KP])+f(i,j,k)*p->DZN[KP]);
            }
        }
    }
    else if(p->A10==4)
    {
        for(n=0;n<gauge_num;++n)
        if(flag[n]>0)
        {
            i = iloc[n];
            j = jloc[n];
            wsf[n] = a->eta(i,j);
        }
    }

    for(n=0;n<gauge_num;++n)
    wsf[n]=pgc->globalmax(wsf[n]);

    // write to file
    if(p->mpirank==0)
    {
        wsfout<<std::setprecision(precision)<<p->simtime<<"\t";
        for(n=0;n<gauge_num;++n)
        {
            wsfout<<std::setprecision(precision)<<wsf[n]<<"\t";
            // flush print to disc limited to prevent data loss for many gauges
            if(n%fileFlushMaxCount==0&&n!=0)
            wsfout<<std::flush;
        }
        wsfout<<std::endl;
    }
}

void print_wsf::ini_location(lexer *p)
{
    for(n=0;n<gauge_num;++n)
    {
        iloc[n] = p->posc_i(x[n]);
        jloc[n] = (p->j_dir ? p->posc_j(y[n]) : 0);

        if(1==ij_boundcheck(p,iloc[n],jloc[n],0))
        flag[n] = 1;
    }
}
