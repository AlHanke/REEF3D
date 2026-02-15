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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void ioflow_f::phi_relax(lexer *p, ghostcell *pgc, field &f)
{
    double relax,distot,distcount;


    if(p->B71>0 && p->count==0)
    LOOP
    {
        distot = 0.0;
        distcount=0;
        for(n=0;n<p->B71;++n)
        {
            dist_B71[n] =  distcalc(p,p->B71_x[n],p->B71_y[n],tan_betaB71[n]);

            if(dist_B71[n]<p->B71_dist[n])
            {
                val = f(i,j,k);
                f(i,j,k)=0.0;
                distot += dist_B71[n];
                ++distcount;
            }
        }


        for(n=0;n<p->B71;++n)
        {
            if(dist_B71[n]<p->B71_dist[n])
            {
                relax = r1(p,dist_B71[n],p->B71_dist[n]);

                if(distcount==1)
                f(i,j,k) += (1.0-relax)*(p->B71_val[n]-p->pos_z()) + relax*val;
                else if(distcount>1)
                f(i,j,k) += ((1.0-relax)*(p->B71_val[n]-p->pos_z()) + relax*val) * (1.0 - dist_B71[n]/(distot>1.0e-10?distot:1.0e20));
            }
        }
    }
}

double ioflow_f::r1(lexer *p, double x, double threshold)
{
    double r=0.0;

    x=1.0-x/(fabs(threshold)>1.0e-10?threshold:1.0e20);
    x=MAX(x,0.0);
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(exp(1.0)-1.0);

    return r;
}

double ioflow_f::distcalc(lexer *p,double x0, double y0, double tan_beta)
{
    double x1,y1;
    double dist=1.0e20;

    x1 = p->pos_x();
    y1 = p->pos_y();

    dist = fabs(y1 - tan_beta*x1 + tan_beta*x0 - y0)/sqrt(pow(tan_beta,2.0)+1.0);

    return dist;
}

int ioflow_f::iozonecheck(lexer *p, fdm*a)
{
    return 1;
}
