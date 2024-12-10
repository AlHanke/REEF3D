/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void partres::timestep(lexer *p, ghostcell *pgc)
{
    double maxVelU=.00,maxVelV=0.0,maxVelW=0.0;
    double maxvz=0.0;
    
    for(size_t n=0;n<P.index;n++)
    {
        if(P.Flag[n]>=0)
        {
            maxVelU = std::max(maxVelU,fabs(P.U[n]));
            maxVelV = std::max(maxVelV,fabs(P.V[n]));
            maxVelW = std::max(maxVelW,fabs(P.W[n]));
        }
    }
    
    
    maxvz = std::max(maxVelU,maxVelV);
    maxvz = std::max(maxvz,maxVelV);
    
    maxvz = pgc->globalmax(maxvz);
    
    maxVelU = pgc->globalmax(maxVelU);
    maxVelV = pgc->globalmax(maxVelV);
    maxVelW = pgc->globalmax(maxVelW);
    
    if(!timestep_ini)
    {
    maxvz = std::max(maxvz,1000.0);
    timestep_ini=true;
    }
    
    
    if(p->S15==0)
        p->dtsed=std::min(p->S13, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));

    if(p->S15==1)
        p->dtsed=std::min(p->dt, (p->S14*p->DXM)/(fabs(maxvz)>1.0e-15?maxvz:1.0e-15));
    
    if(p->S15==2)
    p->dtsed=p->S13;

    p->dtsed=pgc->timesync(p->dtsed);
    p->sedtime+=p->dtsed;
    p->sedtime=pgc->timesync(p->sedtime);
    
    //
	
	if(p->mpirank==0)
    {
	cout<<"tsed: "<<setprecision(4)<<p->sedtime<<" dtsed: "<<setprecision(4)<<p->dtsed<<endl;
    cout<<"Up_max: "<<maxVelU<<endl;
    cout<<"Vp_max: "<<maxVelV<<endl;
    cout<<"Wp_max: "<<maxVelW<<endl;
    }
    
    
}