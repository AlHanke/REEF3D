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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::ray_cast_z(lexer *p, fdm_nhf *d, ghostcell *pgc, int ts, int te)
{
	double ys,ye,zs,ze;
	double Px,Py,Pz;
	double Qx,Qy,Qz;
	double Rx,Ry,Rz;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Cx,Cy,Cz;
	double PQx,PQy,PQz;
	double Mx,My,Mz;
	int is,ie,js,je,ks,ke;
    int checkin;
	double u,v,w;
	double denom;
	double psi = 1.0e-8*p->DXM;

	for(n=ts; n<te; ++n)
	{ 
	Ax = tri_x[n][0];
	Ay = tri_y[n][0];
	Az = tri_z[n][0];
		
	Bx = tri_x[n][1];
	By = tri_y[n][1];
	Bz = tri_z[n][1];
		
	Cx = tri_x[n][2];
	Cy = tri_y[n][2];
	Cz = tri_z[n][2];
	
	
    checkin = 0;
    
	if(Ax>=p->global_xmin && Ax<=p->global_xmax 
    && Ay>=p->global_ymin && Ay<=p->global_ymax
    && Az>=p->global_zmin && Az<=p->global_zmax)
    checkin=1;
    
    if(Bx>=p->global_xmin && Bx<=p->global_xmax 
    && By>=p->global_ymin && By<=p->global_ymax
    && Bz>=p->global_zmin && Bz<=p->global_zmax)
    checkin=1;
    
    if(Cx>=p->global_xmin && Cx<=p->global_xmax 
    && Cy>=p->global_ymin && Cy<=p->global_ymax
    && Cz>=p->global_zmin && Cz<=p->global_zmax)
    checkin=1;
    
    checkin=1;
        
    if(checkin==1)
    {
	xs = MIN3(Ax,Bx,Cx); 
	xe = MAX3(Ax,Bx,Cx);
	
	ys = MIN3(Ay,By,Cy);
	ye = MAX3(Ay,By,Cy);
	
	is = p->posc_i(xs);
	ie = p->posc_i(xe);
	
	js = p->posc_j(ys);
	je = p->posc_j(ye);
		
	
    xs = MIN3(Ax,Bx,Cx) - epsi*p->DXP[is + marge];
	xe = MAX3(Ax,Bx,Cx) + epsi*p->DXP[ie + marge];
	
	ys = MIN3(Ay,By,Cy) - epsi*p->DYP[js + marge];
	ye = MAX3(Ay,By,Cy) + epsi*p->DYP[je + marge];

	
	is = p->posc_i(xs);
	ie = p->posc_i(xe);
	
	js = p->posc_j(ys);
	je = p->posc_j(ye);
	
	is = MAX(is,0);
	ie = MIN(ie,p->knox);
	
	js = MAX(js,0);
	je = MIN(je,p->knoy);
	
		for(i=is;i<ie;i++)
		for(j=js;j<je;j++)
		{
		Px = p->XP[IP]-psi;
		Py = p->YP[JP]+psi;
		Pz = p->global_zmin-10.0*p->DXM ;
		
		Qx = p->XP[IP]+psi;
		Qy = p->YP[JP]-psi;
		Qz = p->global_zmax+10.0*p->DXM ;

		
		PQx = Qx-Px;
		PQy = Qy-Py;
		PQz = Qz-Pz;
		
		// uvw
		Mx = PQy*Pz - PQz*Py;
		My = PQz*Px - PQx*Pz;
		Mz = PQx*Py - PQy*Px;

		
		u = PQx*(Cy*Bz - Cz*By) + PQy*(Cz*Bx - Cx*Bz) + PQz*(Cx*By - Cy*Bx)
		  + Mx*(Cx-Bx) + My*(Cy-By) + Mz*(Cz-Bz);
		  
		v = PQx*(Ay*Cz - Az*Cy) + PQy*(Az*Cx - Ax*Cz) + PQz*(Ax*Cy - Ay*Cx)
		  + Mx*(Ax-Cx) + My*(Ay-Cy) + Mz*(Az-Cz);
		  
		w = PQx*(By*Az - Bz*Ay) + PQy*(Bz*Ax - Bx*Az) + PQz*(Bx*Ay - By*Ax)
		  + Mx*(Bx-Ax) + My*(By-Ay) + Mz*(Bz-Az);
        
        int check=1;
		if(fabs(u)<=1.0e-20 && fabs(v)<=1.0e-20 && fabs(w)<=1.0e-20)
		check = 0;

			if(((u>1.0e-20 && v>1.0e-20 && w>1.0e-20) || (u<-1.0e-20 && v<-1.0e-20 && w<-1.0e-20)) && check==1)
			{
			denom = 1.0/(u+v+w);
            
			u *= denom;
			v *= denom;
			w *= denom;
			
			Rz = u*Az + v*Bz + w*Cz;
            
             k = p->posc_sig(i,j,Rz);

            int distcheck=1;
  
            if(Rz<p->ZSP[IJK])
            if(k>=0 && k<p->knoz)
            if(IO[IJK]<0 && IO[IJKm1]<0)
            distcheck=0;
            
            if(Rz>=p->ZSP[IJK])
            if(k>=0 && k<p->knoz)
            if(IO[IJK]<0 && IO[IJKp1]<0)
            distcheck=0;
            
            if(distcheck==1)
            for(k=0;k<p->knoz;++k)
            d->SOLID[IJK]=MIN(fabs(Rz-p->ZSP[IJK]),d->SOLID[IJK]);
            }
            
		}
	}
    }

}
