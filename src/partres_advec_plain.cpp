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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

    /**
     * @brief Moves the particles with the flow field.
     *
     * This function is responsible for moving the particle with the flow field.
     * It does so by calculating the experienced accelerations and the resulting new velocity and position of the particle.
     *
     * @param p A pointer to the lexer object.
     * @param a A reference to the fdm object.
     * @param pgc A reference to the ghostcell object.
     * @param PP A reference to the particles_obj object.
     */
     
void partres::advec_plain(lexer *p, fdm &a, particles_obj2::particle &_particle, double d50, sediment_fdm &s, turbulence &pturb, 
                        double PX, double PY, double PZ, double PU, double PV, double PW,
                        double &du, double &dv, double &dw, double alpha)
{
    double RKu,RKv,RKw;
    double u=0,v=0,w=0;
    double du1, du2, du3, dv1, dv2, dv3, dw1, dw2, dw3;
    double ws=0;
    double topoDist=0;
       
    double pressureDivX=0, pressureDivY=0, pressureDivZ=0;
    double stressDivX=0, stressDivY=0, stressDivZ=0;
    double netBuoyX=0, netBuoyY=0, netBuoyZ=0;
    // double netBuoyX=(1.0-drho)*p->W20, netBuoyY=(1.0-drho)*p->W21, netBuoyZ=(1.0-drho)*p->W22;

    bool limited = true;
    bool debugPrint = false;
    bool bedLoad = false;
    bool shearVel = true;

    double Urel=0.0,Vrel=0.0,Wrel=0.0;
    double thetas=0;
    double DragCoeff=0;
    double Uabs_rel=0.0;

    i=p->posc_i(PX);
    j=p->posc_j(PY);
    k=p->posc_k(PZ);
        

    velDist=0.6;
    
    u=p->ccipol1(a.u,PX,PY,PZ+velDist*p->DZP[KP]);
    v=p->ccipol2(a.v,PX,PY,PZ+velDist*p->DZP[KP]);


    // relative velocity
    Urel=u-PU;
    Vrel=v-PV;

    _particle.uf=Urel;
    _particle.vf=Urel;
    
    // particle force
    if(p->Q202==1)
    {
    DragCoeff=drag_model(p,d50,Urel,Vrel,Wrel,thetas);
    
    du=DragCoeff*Urel;
    dv=DragCoeff*Vrel;
    }
        
    if(p->Q202==2)
    {
    //relative_velocity(p,a,PP,n,Urel,Vrel,Wrel);
    Uabs_rel=sqrt(Urel*Urel + Vrel*Vrel);
    Re_p = Uabs_rel*d50/(p->W2/p->W1);
    DragCoeff = drag_coefficient(Re_p);
         
    // acceleration
    Fd = p->W1 * DragCoeff * PI/8.0 * pow(d50,2)  * pow(Uabs_rel,2.0);
    
    // relax
    Fd *= rf(p,PX,PY);
    
    // resistance force
    Fs = p->Q30*(p->S22-p->W1)*fabs(p->W22)*PI*pow(d50, 3.0)/6.0;
     
    F_tot = Fd-Fs;//*s.reduce(i,j);
    
    F_tot = MAX(F_tot,0.0);
    
    //cout<<"Fd: "<<Fd<<" Fs: "<<Fs<<" F_tot: "<<F_tot<<" "<<PP.d50<<" "<<DragCoeff<<endl;

    du = F_tot /(p->S22*PI*pow(d50,3.0)/6.0) * Urel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    dv = F_tot /(p->S22*PI*pow(d50,3.0)/6.0) * Vrel/(fabs(Uabs_rel)>1.0e-10?fabs(Uabs_rel):1.0e20);
    
    _particle.shear_eff=Fd;
    _particle.shear_crit=Fs;
    }

    _particle.drag=DragCoeff;
    

    // solid forcing
    double fx,fy,fz;
    /*
    fx = p->ccipol1c(a.fbh1,PX,PY,PZ)*(0.0-PU)/(alpha*p->dtsed); 
    fy = p->ccipol2c(a.fbh2,PX,PY,PZ)*(0.0-PV)/(alpha*p->dtsed); 
    fz = p->ccipol3c(a.fbh3,PX,PY,PZ)*(0.0-PW)/(alpha*p->dtsed); 
    
    du += fx;
    dv += fy;
    dw += fz;*/
    
    // relax
    du *= rf(p,PX,PY);
    dv *= rf(p,PX,PY);
    dw *= rf(p,PX,PY);
    
    dw=0.0;
    
    Umax = MAX(Umax,sqrt(PU*PU + PV*PV));
    
    
    if(PU!=PU || PV!=PV || PW!=PW)
    {
    cout<<"NaN detected.\nUrel: "<<Urel<<" Vrel: "<<Vrel<<" Wrel: "<<Wrel<<"\nDrag: "<<DragCoeff<<endl;
    exit(1);
    }
}
