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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ifou.h"
#include"lexer.h"
#include"fdm.h"

ifou::ifou(lexer *p)
{
    if(p->j_dir==0)
    {
        if(p->B200>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux.emplace<flux_face_FOU_vrans_2D>(p);

            else if(p->D11==2)
            pflux.emplace<flux_face_CDS2_vrans_2D>();
        }
        else
        {
            if(p->D11==1)
            pflux.emplace<flux_face_FOU_2D>(p);

            else if(p->D11==2)
            pflux.emplace<flux_face_CDS2_2D>();
        }
    }
    else if(p->j_dir==1)
    {
        if(p->B200>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux.emplace<flux_face_FOU_vrans>(p);

            else if(p->D11==2)
            pflux.emplace<flux_face_CDS2_vrans>();
        }
        else
        {
            if(p->D11==1)
            pflux.emplace<flux_face_FOU>(p);

            else if(p->D11==2)
            pflux.emplace<flux_face_CDS2>();
        }
    }
}

void ifou::start(lexer *p, fdm *a, field &b, int ipol, field &uvel, field &vvel, field &wvel)
{
    count=0;
    std::visit([&](auto& flux)
    {
        if(ipol==1)
        {
            FIELDLOOP_INC_MEMBER(a,F,
                FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
                aij(flux,p,a,b,1,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZN.data());
            )
        }
        else if(ipol==2 && p->j_dir==1)
        {
            FIELDLOOP_INC_MEMBER(a,G,
                FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
                aij(flux,p,a,b,2,uvel,vvel,wvel,p->DXN.data(),p->DYP.data(),p->DZN.data());
            )
        }
        else if(ipol==3)
        {
            FIELDLOOP_INC_MEMBER(a,H,
                FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
                aij(flux,p,a,b,3,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZP.data());
            )
        }
        else if(ipol==4)
        {
            FIELDLOOP_INC_MEMBER(a,L,
                FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
                aij(flux,p,a,b,4,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZN.data());
            )
        }
    }, pflux);
}

template<typename FluxT, typename GenericField>
void ifou::aij(FluxT &pflux, lexer *p, fdm *a, const GenericField &b, int ipol, const GenericField &uvel, const GenericField &vvel, const GenericField &wvel, double *DX, double *DY, double *DZ)
{
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    pflux.u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux.v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux.w_flux(a,ipol,wvel,kvel1,kvel2);

    double udir,vdir,wdir;
    udir=vdir=wdir=0.0;
    if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;

    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;

    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;

    const double invDXM1=1.0/DX[IM1], invDXP=1.0/DX[IP], invDYM1=1.0/DY[JM1], invDYP=1.0/DY[JP], invDZM1=1.0/DZ[KM1], invDZP=1.0/DZ[KP];

    a->M.p[count] =    udir*ivel2*invDXM1 - (1.0-udir)*ivel1*invDXP
                  + (p->j_dir?(vdir*jvel2*invDYM1 - (1.0-vdir)*jvel1*invDYP):0.0)
                  +  wdir*kvel2*invDZM1 - (1.0-wdir)*kvel1*invDZP;

    a->M.s[count] = -udir*ivel1*invDXM1;
    a->M.n[count] =  (1.0-udir)*ivel2*invDXP;

    if(p->j_dir)
    {
        a->M.e[count] = -vdir*jvel1*invDYM1;
        a->M.w[count] = (1.0-vdir)*jvel2*invDYP;
    }
    else
    {
        a->M.e[count] = 0.0;
        a->M.w[count] = 0.0;
    }

    a->M.b[count] = -wdir*kvel1*invDZM1;
    a->M.t[count] =  (1.0-wdir)*kvel2*invDZP;

    ++count;
}
