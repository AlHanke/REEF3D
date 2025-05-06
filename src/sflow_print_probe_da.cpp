/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"sflow_print_probe_da.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<string>
#include<iostream>

sflow_print_probe_da::sflow_print_probe_da(lexer *p, fdm2D *b, ghostcell *pgc) : probenum(p->P63)
{
    p->Iarray(iloc,probenum);
    p->Iarray(jloc,probenum);
    p->Iarray(flag,probenum);
    
    if(p->mpirank==0 && probenum>0)
    {
        cout<<"probepoint_num: "<<probenum<<endl;

        // Create Folder
        mkdir("./REEF3D_SFLOW_ProbePoint",0777);

        pout = new ofstream[probenum];
        // write file headers
        for(int n=0;n<probenum;++n)
        {
            char name[100];
            sprintf(name,"./REEF3D_SFLOW_ProbePoint/REEF3D-SFLOW-Probe-Point-%i.dat",n+1);
            
            pout[n].open(name);

            pout[n]<<"Depth Averaged Point Probe ID:  "<<n<<"\n\n";
            pout[n]<<"x_coord     y_coord\n";
            
            pout[n]<<n+1<<"\t "<<p->P63_x[n]<<"\t "<<p->P63_y[n]<<"\n\n\n";
            
            pout[n]<<"t \t Um \t Vm \t Wm \t Pm \t eta\n";

            pout[n].close();
        }
    }
    
    ini_location(p);
}

sflow_print_probe_da::~sflow_print_probe_da()
{
    delete[] iloc;
    delete[] jloc;
    delete[] flag;
    
    delete[] pout;
}

void sflow_print_probe_da::start(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // find values and write
    for(int n=0;n<probenum;++n)
    {
        double uval = -1.0e20;
        double vval = -1.0e20;
        double wval = -1.0e20;
        double pval = -1.0e20;
        double eval = -1.0e20;
    
        if(flag[n]>0)
        {
            double xp =p->P63_x[n];
            double yp = p->P63_y[n];
            
            uval = p->ccslipol1(b->P, xp, yp);
            vval = p->ccslipol2(b->Q, xp, yp);
            wval = p->ccslipol4(b->ws,xp,yp);
            pval = p->ccslipol4(b->press,xp,yp);
            eval = p->ccslipol4(b->eta,xp,yp);
        }
    
        uval = pgc->globalmax(uval);
        vval = pgc->globalmax(vval);
        wval = pgc->globalmax(wval);
        pval = pgc->globalmax(pval);
        eval = pgc->globalmax(eval);
    
        if(p->mpirank==0)
        {
            std::string name = "./REEF3D_SFLOW_ProbePoint/REEF3D-SFLOW-Probe-Point-" + std::to_string(n + 1) + ".dat";
            pout[n].open(name, std::ofstream::out | std::ofstream::app);
            pout[n]<<setprecision(9)<<p->simtime<<" \t "<<uval<<" \t "<<vval<<" \t "<<wval<<" \t "<<pval<<" \t "<<eval<<"\n";
            pout[n].close();
        }
    }
}

void sflow_print_probe_da::ini_location(lexer *p)
{
    for(int n=0;n<probenum;++n)
    {
        bool check=false;
        
        iloc[n]=p->posc_i(p->P63_x[n]);
        
        if(p->j_dir==0)
            jloc[n]=0;
        else if(p->j_dir==1)
            jloc[n]=p->posc_j(p->P63_y[n]);
        
        if(iloc[n]>=0 && iloc[n]<p->knox)
            if(jloc[n]>=0 && jloc[n]<p->knoy)
                check=true;
        
        if(check)
        {
            int i = iloc[n];
            int j = jloc[n];
            
            if(p->flagslice4[IJ]<0)
                check=false;
        }

        if(check)
            flag[n]=1;
    }
}
