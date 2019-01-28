/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"fnpf_sg_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"field4.h"
#include"convection.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"onephase.h"
#include"fnpf_voiddisc.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"fnpf_cds6.h"
#include"fnpf_weno.h"
#include"fnpf_weno_wd.h"
#include"fnpf_wenoflux.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"

fnpf_sg_fsfbc::fnpf_sg_fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc) : diss(p), Fxx(p), Fyy(p)
{    
    if(p->A311==0)
    pconvec = new fnpf_voiddisc(p);
    
    if(p->A311==2)
    pconvec = new fnpf_cds2(p);
    
    if(p->A311==3)
    pconvec = new fnpf_cds4(p);
    
    if(p->A311==4)
    pconvec = new fnpf_weno(p);
    
    if(p->A311==5)
    pconvec = new fnpf_weno_wd(p,c);
    
    if(p->A311==6)
    pconvec = new fnpf_cds6(p);
    
    pdh = new fnpf_weno(p);
    
    if(p->A312==2)
    {
    pddx = new fnpf_ddx_cds2(p);
    pdx = new fnpf_cds2(p);
    }
    
    if(p->A312==3)
    {
    pddx = new fnpf_ddx_cds4(p);
    pdx = new fnpf_cds4(p);
    }
    
    
    FFILOOP4
    {
    c->Fy(i,j) = 0.0;
    c->Ey(i,j) = 0.0;
    c->Hy(i,j) = 0.0;
    c->Eyy(i,j) = 0.0;
    Fyy(i,j) = 0.0;
    }
    
    
    wd_criterion=0.00005;
    
    if(p->A344==1)
    wd_criterion=p->A244_val;
    
    if(p->A345==1)
    wd_criterion=p->A245_val*p->dx;
}

fnpf_sg_fsfbc::~fnpf_sg_fsfbc()
{
}

void fnpf_sg_fsfbc::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = MAX(0.0,c->eta(i,j) + p->wd - c->bed(i,j));
    
    // fi
    if(p->i_dir==1 && p->j_dir==1)
    FFILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
    jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
    c->Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    c->Fy(i,j) = pconvec->sy(p,Fifsf,jvel);
    
    c->Ex(i,j) = pdh->sx(p,eta,ivel);
    c->Ey(i,j) = pdh->sy(p,eta,jvel);
    
    c->Exx(i,j) = pddx->sxx(p,eta);
    c->Eyy(i,j) = pddx->syy(p,eta);
    
    Fxx(i,j) = pddx->sxx(p,Fifsf);
    Fyy(i,j) = pddx->syy(p,Fifsf);
    }
    
    if(p->i_dir==1 && p->j_dir==0)
    FFILOOP4
    {
    ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);    
    
    c->Fx(i,j) = pconvec->sx(p,Fifsf,ivel);
    c->Ex(i,j) = pdh->sx(p,eta,ivel);
    c->Exx(i,j) = pddx->sxx(p,eta);
    Fxx(i,j) = pddx->sxx(p,Fifsf);
    }

}

void fnpf_sg_fsfbc::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    FFILOOP4
    {
    c->Bx(i,j) = pdx->sx(p,c->depth,1.0);
    c->By(i,j) = pdx->sy(p,c->depth,1.0);
    
    c->Bxx(i,j) = pddx->sxx(p,c->depth);
    c->Byy(i,j) = pddx->syy(p,c->depth);
    }
    
    pgc->gcsl_start4(p,c->Bx,1);
    pgc->gcsl_start4(p,c->By,1);
    
    SLICELOOP4
    c->wet(i,j)=1;
    
    pgc->gcsl_start4int(p,c->wet,50);
}

void fnpf_sg_fsfbc::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    FFILOOP4
    c->Fz(i,j) = p->sigz[IJ]*pconvec->sz(p,c->Fi);
}

void fnpf_sg_fsfbc::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
    {
    c->K(i,j) =  - c->Fx(i,j)*c->Ex(i,j) - c->Fy(i,j)*c->Ey(i,j) 
    
                 + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0)) 
                 
                 + diss(i,j)*(c->Exx(i,j) + c->Eyy(i,j));
    }
}

void fnpf_sg_fsfbc::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{
    SLICELOOP4
    c->K(i,j) =  - 0.5*c->Fx(i,j)*c->Fx(i,j) - 0.5*c->Fy(i,j)*c->Fy(i,j) 
    
                 + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j)
                 
                 + diss(i,j)*(Fxx(i,j) + Fyy(i,j));
}


