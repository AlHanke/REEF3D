/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sigma.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_sg_fsfbc.h"

fnpf_sigma::fnpf_sigma(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
}

fnpf_sigma::~fnpf_sigma()
{
}

void fnpf_sigma::sigma_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_sg_fsfbc *pf, slice &eta)
{	
    p->Darray(p->sig,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigz,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+1));
    
    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sig[FIJKm1] = p->ZN[KM1];
            p->sig[FIJKm2] = p->ZN[KM2];
            p->sig[FIJKm3] = p->ZN[KM3];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sig[FIJKp1] = p->ZN[KP1];
            p->sig[FIJKp2] = p->ZN[KP2];
            p->sig[FIJKp3] = p->ZN[KP3];
        } 
    }
    
    
    SLICELOOP4
	c->bed(i,j) = p->bed[IJ];
    
    SLICELOOP4
    c->WL(i,j) = c->eta(i,j) + p->wd - c->bed(i,j);
    
    SLICELOOP4
	c->depth(i,j) = p->wd - c->bed(i,j);
    
    pgc->gcsl_start4(p,c->depth,50);
}

void fnpf_sigma::sigma_update(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_sg_fsfbc *pf, slice &eta)
{
    FLOOP
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(pf->Bx(i,j)/c->WL(i,j)) - p->sig[FIJK]*(pf->Ex(i,j)/c->WL(i,j));
    
    FLOOP
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(pf->By(i,j)/c->WL(i,j)) - p->sig[FIJK]*(pf->Ey(i,j)/c->WL(i,j));
    
    FLOOP
    p->sigz[FIJK] = 1.0/c->WL(i,j);
    
    FLOOP
    p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/c->WL(i,j))*(pf->Bxx(i,j) - pow(pf->Bx(i,j),2.0)/c->WL(i,j)) // xx
    
                  - (p->sig[FIJK]/c->WL(i,j))*(pf->Exx(i,j) - pow(pf->Ex(i,j),2.0)/c->WL(i,j))
                  
                  - (p->sigx[FIJK]/c->WL(i,j))*(pf->Bx(i,j) + pf->Ex(i,j))*0.0
                  
                  - ((1.0 - 2.0*p->sig[FIJK])/pow(c->WL(i,j),2.0))*(pf->Bx(i,j)*pf->Ex(i,j))
                  
                  
                  + ((1.0 - p->sig[FIJK])/c->WL(i,j))*(pf->Byy(i,j) - pow(pf->By(i,j),2.0)/c->WL(i,j)) // yy
    
                  - (p->sig[FIJK]/c->WL(i,j))*(pf->Eyy(i,j) - pow(pf->Ey(i,j),2.0)/c->WL(i,j))
                  
                  - (p->sigy[FIJK]/c->WL(i,j))*(pf->By(i,j) + pf->Ey(i,j))*0.0
                  
                  - ((1.0 - 2.0*p->sig[FIJK])/pow(c->WL(i,j),2.0))*(pf->By(i,j)*pf->Ey(i,j));
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigx[FIJKm1] = p->sigx[KP];
            p->sigx[FIJKm2] = p->sigx[KP];
            p->sigx[FIJKm3] = p->sigx[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigx[FIJKp1] = p->sigx[KP];
            p->sigx[FIJKp2] = p->sigx[KP];
            p->sigx[FIJKp3] = p->sigx[KP];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigy[FIJKm1] = p->sigy[KP];
            p->sigy[FIJKm2] = p->sigy[KP];
            p->sigy[FIJKm3] = p->sigy[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigy[FIJKp1] = p->sigy[KP];
            p->sigy[FIJKp2] = p->sigy[KP];
            p->sigy[FIJKp3] = p->sigy[KP];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigxx[FIJKm1] = p->sigxx[KP];
            p->sigxx[FIJKm2] = p->sigxx[KP];
            p->sigxx[FIJKm3] = p->sigxx[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigxx[FIJKp1] = p->sigxx[KP];
            p->sigxx[FIJKp2] = p->sigxx[KP];
            p->sigxx[FIJKp3] = p->sigxx[KP];
        } 
    }
    
    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*(c->eta(i,j) + p->wd - c->bed(i,j)) + c->bed(i,j);
    

}




