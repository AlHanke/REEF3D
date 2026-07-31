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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"weno3_flux.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"

weno3_flux::weno3_flux(lexer* p) : weno3_nug_func(p)
{
    if(p->j_dir==0)
    {
        if(p->B200==0)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_2D(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_2D;
        }
        else if(p->B200>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans_2D(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans_2D;
        }
    }
    else if(p->j_dir==1)
    {
        if(p->B200==0)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2;
        }
        else if(p->B200>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans;
        }
    }
}

void weno3_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;

    if(ipol==1)
    {
        uf=1;
        FIELDLOOP_INC_MEMBER(
            a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZN.data());
        )
    }
    else if(ipol==2 && p->j_dir==1)
    {
        vf=1;
        FIELDLOOP_INC_MEMBER(
            a,G,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN.data(),p->DYP.data(),p->DZN.data());
        )
    }
    else if(ipol==3)
    {
        wf=1;
        FIELDLOOP_INC_MEMBER(
            a,H,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZP.data());
        )
    }
    else if(ipol==4)
    {
        FIELDLOOP_INC_MEMBER(
            a,L,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            L(i,j,k) += aij(p,a,b,4,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZN.data());
        )
    }
}

template<typename GenericField>
double weno3_flux::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel, double *DX, double *DY, double *DZ)
{
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    const double du = fx_div(b,i,j,k,ivel1,ivel2,IP,uf,epsilon,psi);
    const double dw = fz_div(b,i,j,k,kvel1,kvel2,KP,wf,epsilon,psi);

    if(p->j_dir==1)
    {
        const double dv = fy_div(b,i,j,k,jvel1,jvel2,JP,vf,epsilon,psi);
        return -(du/DX[IP] + dv/DY[JP] + dw/DZ[KP]);
    }
    return -(du/DX[IP] + dw/DZ[KP]);
}
