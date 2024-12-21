/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"wave_interface.h"
#include"wave_lib_void.h"
#include"wave_lib_shallow.h"
#include"wave_lib_deep.h"
#include"wave_lib_linear.h"
#include"wave_lib_flap.h"
#include"wave_lib_flap_double.h"
#include"wave_lib_piston.h"
#include"wave_lib_piston_eta.h"
#include"wave_lib_flap_eta.h"
#include"wave_lib_Stokes_2nd.h"
#include"wave_lib_Stokes_5th.h"
#include"wave_lib_Stokes_5th_SH.h"
#include"wave_lib_cnoidal_shallow.h"
#include"wave_lib_cnoidal_1st.h"
#include"wave_lib_cnoidal_5th.h"
#include"wave_lib_solitary_1st.h"
#include"wave_lib_solitary_3rd.h"
#include"wave_lib_irregular_1st.h"
#include"wave_lib_irregular_2nd_a.h"
#include"wave_lib_irregular_2nd_b.h"
#include"wave_lib_reconstruct.h"
#include"wave_lib_hdc.h"
#include"wave_lib_ssgw.h"
#include"lexer.h"

wave_interface::wave_interface(lexer *p, ghostcell *pgc) 
{ 
    p->wts=0.0;
    p->wte=1.0e20;
    
    if(p->B94==0)
	 wD=p->phimean;
	
	if(p->B94==1)
	wD=p->B94_wdt;

    switch(p->B92) // wave type
    {
    default:
    case 0:
        pwave = new wave_lib_void(p);
        break;
	
    case 1:
        pwave = new wave_lib_shallow(p);
        break;
        
    case 2:
        pwave = new wave_lib_linear(p);
        break;
        
    case 3:
        pwave = new wave_lib_deep(p);
        break;
        
    case 4:
        pwave = new wave_lib_Stokes_2nd(p);
        break;
        
    case 5:
        pwave = new wave_lib_Stokes_5th(p);
        break;
        
    case 6:
        pwave = new wave_lib_cnoidal_shallow(p);
        break;
        
    case 7:
        pwave = new wave_lib_cnoidal_1st(p);
        break;
        
    case 8:
        pwave = new wave_lib_cnoidal_5th(p);
        break;
        
    case 9:
        pwave = new wave_lib_solitary_1st(p);
        break;
        
    case 10:
        pwave = new wave_lib_solitary_3rd(p);
        break;
        
    case 11:
        pwave = new wave_lib_Stokes_5th_SH(p);
        break;
        
    case 20:
        pwave = new wave_lib_piston_eta(p);
        break;
        
    case 21:
        pwave = new wave_lib_piston(p);
        break;
        
    case 22:
        pwave = new wave_lib_flap(p);
        break;
        
    case 23:
        pwave = new wave_lib_flap_double(p);
        break;
        
    case 24:
        pwave = new wave_lib_flap_eta(p);
        break;
        
    case 31:
    case 41:
    case 51:
        pwave = new wave_lib_irregular_1st(p,pgc);
        break;
        
    case 32:
    case 42:
    case 52:
        pwave = new wave_lib_irregular_2nd_a(p,pgc);
        break;
        
    case 33:
    case 43:
    case 53:
        pwave = new wave_lib_irregular_2nd_b(p,pgc);
        break;
        
    case 61:
        pwave = new wave_lib_hdc(p);
        break;
    
    case 70:
        pwave = new wave_lib_ssgw(p);
        break;
    }
}

double wave_interface::wave_u(lexer *p, ghostcell *, double x, double y, double z)
{

    double uvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u(p,x,y,z);
	
    return uvel;
}

double wave_interface::wave_v(lexer *p, ghostcell *, double x, double y, double z)
{
    double vvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    vvel = pwave->wave_v(p,x,y,z);

    return vvel;
}

double wave_interface::wave_w(lexer *p, ghostcell *, double x, double y, double z)
{
    double wvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    wvel = pwave->wave_w(p,x,y,z);

    return wvel;
}

double wave_interface::wave_h(lexer *p, ghostcell *, double x, double y, double )
{
    double lsv=p->phimean;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    lsv=p->phimean + pwave->wave_eta(p,x,y);

    return lsv;
}

double wave_interface::wave_fi(lexer *p, ghostcell *, double x, double y, double z)
{
    double pval=0.0;
    
    z = MAX(z,-wD);
    
    pval = pwave->wave_fi(p,x,y,z);
	
    return pval;
}

double wave_interface::wave_eta(lexer *p, ghostcell *, double x, double y)
{
    double eta=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    eta = pwave->wave_eta(p,x,y);
	
    return eta;
}

void wave_interface::wave_prestep(lexer *p)
{
    pwave->wave_prestep(p);
}

int wave_interface::printcheck=0;
