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
#include "wsf_locate.h"
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

    lev.resize(gauge_num);
    levmax.resize(gauge_num);
    wsf.resize(gauge_num);
}

print_wsf::~print_wsf()
{
    wsfout.close();
}

void print_wsf::height_gauge(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    std::fill(wsf.begin(), wsf.end(), -1.0e20);
    std::fill(lev.begin(), lev.end(), -1);

    // Gauge locations are resolved here, not at construction: printer_CFD is
    // built in logic_cfd(), which runs before driver_ini_cfd()'s regrid loop,
    // so a location fixed at construction only ever saw a single-level grid.
    // Re-resolving every call also tracks a surface that moves between levels
    // as the hierarchy is re-tagged during the run.
    LEVEL_LOOP
    TILE_LOOP
    {
        for(n=0;n<gauge_num;++n)
        {
            if(!locate(p,n,i,j))
            continue;

            if(p->A10==6 && p->F80!=4)
            {
                KLOOP
                PCHECK
                {
                    if(f(i,j,k)>=0.0 && f(i,j,k+1)<0.0)
                    record(n,-(f(i,j,k)*p->DZP[KP])/(f(i,j,k+1)-f(i,j,k)) + p->pos_z());
                }
            }
            else if(p->A10==6 && p->F80==4)
            {
                KLOOP
                {
                    if(f(i,j,k)>p->F94 && f(i,j,k+1)<p->F93)
                    record(n,p->pos_z()+0.5*p->DZN[KP]);

                    else if(f(i,j,k)<=p->F94 && f(i,j,k)>=p->F93)
                    record(n,(p->pos_z()-0.5*p->DZN[KP])+f(i,j,k)*p->DZN[KP]);
                }
            }
            else if(p->A10==4)
            record(n,a->eta(i,j));
        }
    }

    // Two collectives total, not two per gauge: agree elementwise on the finest
    // level that produced a surface anywhere, then blank every contribution from
    // a coarser one so the value reduction cannot pick a coarse answer over a
    // fine one. levmax is a copy because the reduction is in place and the local
    // level is still needed to test against the reduced one. Where no rank
    // contributed, levmax stays -1, nothing is blanked, and the gauge reports the
    // sentinel as before.
    levmax = lev;
    pgc->globalimax(levmax.data(),gauge_num);

    for(n=0;n<gauge_num;++n)
    if(lev[n]!=levmax[n])
    wsf[n] = -1.0e20;

    pgc->globalmax(wsf.data(),gauge_num);

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

void print_wsf::record(int gauge, double value)
{
    if(level>lev[gauge])
    {
        lev[gauge] = level;
        wsf[gauge] = value;
    }

    else if(level==lev[gauge])
    wsf[gauge] = MAX(wsf[gauge],value);
}

bool print_wsf::locate(lexer *p, int gauge, int& ii, int& jj) const
{
    ii = wsf_locate::cell_1d(p->XN,ORIGIN_I,IMAX_LOOP,x[gauge]);

    if(ii<0)
    return false;

    if(p->j_dir)
    {
        jj = wsf_locate::cell_1d(p->YN,ORIGIN_J,JMAX_LOOP,y[gauge]);

        if(jj<0)
        return false;
    }

    else
    jj = 0;

    return true;
}
