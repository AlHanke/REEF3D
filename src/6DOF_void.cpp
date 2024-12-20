/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"6DOF_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include<sys/stat.h>

#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"
#include"net.h"
#include"net_void.h"
#include"net_barDyn.h"
#include"net_barQuasiStatic.h"
#include"net_sheet.h"
    
sixdof_void::sixdof_void(lexer* p)
{
    if(p->mpirank==0)
    mkdir("./REEF3D_CFD_6DOF",0777);
}

void sixdof_void::ini(lexer*, ghostcell*)
{
}

void sixdof_void::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int , field &, field &, field &, field &, field &, field &, bool )
{

    if (p->X310 > 0)
    {
        for (int i=0; i<p->mooring_count; i++)
        {
            pmooring[i]->start(p, pgc);
        }
    }
    
    if (p->X320 > 0)
    {
        for (int ii = 0; ii < p->net_count; ii++)
        {
            pnet[ii]->start(p, a, pgc, 1.0, quatRotMat);
            pvrans->start(p, a, pgc, pnet[ii], ii);

            // Forces on rigid body
            pnet[ii]->netForces(p,Xne[ii],Yne[ii],Zne[ii],Kne[ii],Mne[ii],Nne[ii]);

            if( p->mpirank == 0)
            {
                cout<<"Xne"<< ii <<" : "<<Xne[ii]<<" Yne"<< ii <<" : "<<Yne[ii]<<" Zne"<< ii <<" : "<<Zne[ii]
                <<" Kne"<< ii <<" : "<<Kne[ii]<<" Mne"<< ii <<" : "<<Mne[ii]<<" Nne"<< ii <<" : "<<Nne[ii]<<endl;		
            }
        }
    }
    
    ++p->printcount_sixdof;
}


void sixdof_void::start_nhflow(lexer*, fdm_nhf*, ghostcell* , vrans* , vector<net*>& , int , 
                                        double *, double *, double *, slice &, slice &, bool )
{
}

void sixdof_void::start_sflow(lexer*, ghostcell*, int , slice&, slice &, slice &, slice &, slice &, slice &, slice &, bool )
{
    
}

void sixdof_void::isource(lexer*, fdm*, ghostcell*)
{
}

void sixdof_void::jsource(lexer*, fdm*, ghostcell*)
{
}

void sixdof_void::ksource(lexer*, fdm*, ghostcell*)
{
}

void sixdof_void::isource(lexer*, fdm_nhf*, ghostcell*, slice &)
{
}

void sixdof_void::jsource(lexer*, fdm_nhf*, ghostcell*, slice &)
{
}

void sixdof_void::ksource(lexer*, fdm_nhf*, ghostcell*, slice &)
{
}

void sixdof_void::isource2D(lexer*, fdm2D*, ghostcell*)
{
}

void sixdof_void::jsource2D(lexer*, fdm2D*, ghostcell*)
{
}
