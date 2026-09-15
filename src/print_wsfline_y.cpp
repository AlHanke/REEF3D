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

#include<iomanip>
#include"print_wsfline_y.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"wsf_locate.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<algorithm>

print_wsfline_y::print_wsfline_y(lexer *p, fdm* a, ghostcell *pgc)
{
    yloc.resize(p->P56);
    wsf.resize(p->P56);
    yloc_all.resize(p->P56);
    wsf_all.resize(p->P56);
    wsfpoints.resize(p->P56,0);

    recvcount.resize(p->mpi_size);
    recvdispl.resize(p->mpi_size);

    // Create Folder
    if(p->mpirank==0)
    mkdir("./REEF3D_CFD_WSFLINE_Y",0777);
}

print_wsfline_y::~print_wsfline_y()
{
    wsfout.close();
}

void print_wsfline_y::wsfline(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
    char name[250];
    int num;

    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		sprintf(name,"./REEF3D_CFD_WSFLINE_Y/REEF3D-CFD-wsfline_y-%08i.dat",num);

		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines_y:  "<<p->P56<<endl<<endl;
		wsfout<<"line_No     x_coord"<<endl;
		for(q=0;q<p->P56;++q)
		wsfout<<q+1<<"\t "<<p->P56_x[q]<<endl;

		if(p->P53==1)
		wsfout<<q+1<<"\t "<<" Wave Theory "<<endl;

		wsfout<<endl<<endl;

		wsfout<<"y_coord";
		for(q=0;q<p->P56+p->P53;++q)
		wsfout<<"\t \t P "<<q+1;

		wsfout<<endl<<endl;
    }

    //-------------------

    collect(p,a,pgc);

    for(q=0;q<p->P56;++q)
    assemble(p,pgc,q);

    // write to file
    if(p->mpirank==0)
    {
        // Lines no longer share a point count -- a line crossing a refined patch
        // carries more points than one that does not -- so the row count is the
        // longest line and shorter lines are padded. The old rowflag pass is gone
        // with the fixed-width layout it policed: every gathered record is a real
        // point now, so there is nothing to test for emptiness.
        int maxpoints=0;

        for(q=0;q<p->P56;++q)
        maxpoints = MAX(maxpoints,wsfpoints[q]);

        for(n=0;n<maxpoints;++n)
        {
		    for(q=0;q<p->P56;++q)
			{
				if(n<wsfpoints[q])
				{
				wsfout<<setprecision(5)<<yloc_all[q][n]<<" \t ";
				wsfout<<setprecision(5)<<wsf_all[q][n]<<" \t  ";

					if(p->P53==1)
					wsfout<<pflow->wave_fsf(p,pgc,yloc_all[q][n])<<" \t  ";
				}

				else
				{
				wsfout<<setprecision(5)<<" \t \t  ";
				wsfout<<setprecision(5)<<" \t \t  ";

					if(p->P53==1)
					wsfout<<" \t  ";
				}
			}

            wsfout<<endl;
        }

    wsfout.close();
    }
}

void print_wsfline_y::collect(lexer *p, fdm *a, ghostcell *pgc)
{
    for(q=0;q<p->P56;++q)
    {
        yloc[q].clear();
        wsf[q].clear();
    }

    LEVEL_LOOP TILE_LOOP
    {
        for(q=0;q<p->P56;++q)
        {
            if(!locate_i(p,q,i))
            continue;

            JLOOP
            {
                double zval=-1.0e20;

                KLOOP
                PCHECK
                {
                    if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                    {
                        // Drop the crossing when a finer level overlays this
                        // cell: the coarse centre does not coincide with any
                        // fine centre, so it would survive the dedupe below and
                        // interleave with the fine points instead of being
                        // replaced by them.
                        if(!wsf_locate::uncovered(p,i,j,k))
                        continue;

                        zval=MAX(zval,-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
                    }
                }

                if(zval>-1.0e20)
                {
                    yloc[q].push_back(p->pos_y());
                    wsf[q].push_back(zval);
                }
            }
        }
    }
}

void print_wsfline_y::assemble(lexer *p, ghostcell *pgc, int line)
{
    // Every rank contributes a different number of points, and the number
    // changes with every regrid, so the counts are exchanged each call and the
    // gather is a gatherv. allgather rather than gather: the displacements cost
    // nothing to compute everywhere and the alternative is a second collective.
    int sendcount = int(yloc[line].size());

    pgc->allgather_int(&sendcount,1,recvcount.data(),1);

    int total=0;

    for(int r=0;r<p->mpi_size;++r)
    {
        recvdispl[r]=total;
        total+=recvcount[r];
    }

    if(p->mpirank==0)
    {
        yloc_all[line].assign(size_t(total),0.0);
        wsf_all[line].assign(size_t(total),0.0);
    }

    pgc->gatherv_double(yloc[line].data(),sendcount,yloc_all[line].data(),recvcount.data(),recvdispl.data());
    pgc->gatherv_double(wsf[line].data(),sendcount,wsf_all[line].data(),recvcount.data(),recvdispl.data());

    wsfpoints[line]=total;

    if(p->mpirank==0 && total>0)
    {
        sort(yloc_all[line].data(), wsf_all[line].data(), 0, total-1);
        remove_multientry(p,yloc_all[line].data(), wsf_all[line].data(), wsfpoints[line]);
    }

    if(p->mpirank!=0)
    wsfpoints[line]=0;
}

bool print_wsfline_y::locate_i(lexer *p, int line, int& ii) const
{
    ii = wsf_locate::cell_1d(p->XN,ORIGIN_I,IMAX_LOOP,p->P56_x[line]);

    return ii>=0;
}

void print_wsfline_y::sort(double *a, double *b, int left, int right)
 {

  if (left < right)
  {

    double pivot = a[right];
    int l = left;
    int r = right;

    do {
      while (a[l] < pivot) l++;

      while (a[r] > pivot) r--;

      if (l <= r) {
          double swap = a[l];
          double swapd = b[l];

          a[l] = a[r];
          a[r] = swap;

          b[l] = b[r];
          b[r] = swapd;

          l++;
          r--;
      }
    } while (l <= r);

    sort(a,b, left, r);
    sort(a,b, l, right);
  }
}

void print_wsfline_y::remove_multientry(lexer *p, double* b, double* c, int& num)
{
    int oldnum=num;
    double yval=-1.12e23;

    int count=0;

    std::vector<double> f(size_t(num),0.0);
    std::vector<double> g(size_t(num),-1.12e22);

    // DXM, not DYM, is deliberate and pre-existing: grid::gridspacing overwrites
    // DXM with an all-directions average, so it is the neutral length scale here.
    // The tolerance is a fraction of the COARSE spacing, so it still only merges
    // genuine duplicates (the same column reported by two ranks sharing a halo)
    // and never two distinct fine cells, whose centres are DXM/(2*ref_ratio^lev)
    // apart.
    for(n=0;n<oldnum;++n)
    {
        if(yval<=b[n]+0.001*p->DXM && yval>=b[n]-0.001*p->DXM && count>0)
        g[count-1]=MAX(g[count-1],c[n]);

        if(yval>b[n]+0.001*p->DXM || yval<b[n]-0.001*p->DXM)
        {
        f[count]=b[n];
        g[count]=c[n];
        ++count;
        }

    yval=b[n];
    }

    for(n=0;n<count;++n)
    {
    b[n]=f[n];
    c[n]=g[n];
    }

	num=count;
}
