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

#include"force.h"
#include"lexer.h"
#include<sys/stat.h>
#include<sys/types.h>

void force::print_force(lexer* p)
{
    // write to surf file
    if(sqrt(Fx*Fx+Fy*Fy+Fz*Fz)>=threshold)
    {
        fout<<p->count<<"\t";
        fout<<setprecision(9)<<p->simtime<<"\t";
        fout<<Fx<<" \t "<<Fy<<" \t "<<Fz<<endl;
    }
}

void force::print_ini(lexer* p)
{   
    if(p->mpirank==0)
    {
        // Create Folder
        mkdir("./REEF3D_CFD_Force",0777);

        // open force surf file
        sprintf(name,"./REEF3D_CFD_Force/REEF3D_CFD_Force-%i.dat",ID+1);
        fout.open(name);

        fout<<"x_start xend | y_start y_end | z_start z_end\n";
        
        fout<<p->P81_xs[ID]<<" "<<p->P81_xe[ID]<<" | "<<p->P81_ys[ID]<<" "<<p->P81_ye[ID]<<" | "<<p->P81_zs[ID]<<" "<<p->P81_ze[ID]<<"\n\n\n";

        fout<<"it \t time \t Fx \t Fy \t Fz"<<endl;
    }
}
