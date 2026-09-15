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

#include "print_wsf_theory.h"
#include "lexer.h"
#include "ghostcell.h"
#include "ioflow.h"

#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

print_wsf_theory::print_wsf_theory(lexer *p, fdm*, ghostcell *pgc, int num)
{
    if(p->P50>0 && num==0)
    {
        gauge_num = p->P50;
        x = p->P50_x;
        y = p->P50_y;
    }
    else
    {
        if(p->mpirank==0)
        std::cerr<<"Error: Impropper height gauge defined!"<<std::endl;

        pgc->final(EXIT_FAILURE);
    }

    if(p->mpirank==0 && p->P50>0 && num==0)
    {
        // Create Folder
        mkdir("./REEF3D_CFD_WSF",0777);

        // open file
        wsfout.open("./REEF3D_CFD_WSF/REEF3D-CFD-WSF-HG-THEORY.dat");

        wsfout<<"number of gauges:  "<<gauge_num<<"\n\n";
        wsfout<<"x_coord     y_coord\n";
        for(int n=0; n<gauge_num; ++n)
        wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<"\n\n\n";

        wsfout<<"time";
        for(int n=0; n<gauge_num; ++n)
        wsfout<<"\t P"<<n+1;

        wsfout<<"\n\n"<<std::flush;
    }
}

print_wsf_theory::~print_wsf_theory()
{
    wsfout.close();
}

void print_wsf_theory::height_gauge(lexer *p, fdm*, ghostcell *pgc, ioflow *pflow, field&)
{
    // write to file
    if(p->mpirank==0)
    {
        wsfout<<std::setprecision(9)<<p->simtime<<"\t";
        for(int n=0; n<gauge_num; ++n)
        {
            wsfout<<std::setprecision(9)<<pflow->wave_fsf(p,pgc,x[n]);
            if(n != gauge_num-1)
            wsfout<<"\t";
        }
        wsfout<<std::endl;
    }
}
