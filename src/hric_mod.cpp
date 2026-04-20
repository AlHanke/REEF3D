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

#include"hric_mod.h"
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

hric_mod::hric_mod(lexer *p)
{
    if(p->j_dir==0)
    {
        if(p->B269>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans_2D(p);

            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans_2D;
        }
        else
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_2D(p);

            else if(p->D11==2)
            pflux = new flux_face_CDS2_2D;
        }
    }
    else if(p->j_dir==1)
    {
        if(p->B269>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans(p);

            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans;
        }
        else
        {
            if(p->D11==1)
            pflux = new flux_face_FOU(p);

            else if(p->D11==2)
            pflux = new flux_face_CDS2;
        }
    }
}

inline void hric_mod::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    {
        FIELDLOOP_INC_MEMBER(a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel);
        )
    }
    else if(ipol==2 && p->j_dir==1)
    {
        FIELDLOOP_INC_MEMBER(a,G,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel);
        )
    }
    else if(ipol==3)
    {
        FIELDLOOP_INC_MEMBER(a,H,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel);
        )
    }
    else if(ipol==4)
    {
        FIELDLOOP_INC_MEMBER(a,L,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel);
        )
    }
}

template<typename GenericField>
inline double hric_mod::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    const double fx1 = cface(p,a,b,1,-1,ivel1);
    const double fx2 = cface(p,a,b,1,0,ivel2);
    const double dx = (ivel2*fx2 - ivel1*fx1)/(p->DXM);

    double dy = 0.0;
    if(p->j_dir==1)
    {
        const double fy1 = cface(p,a,b,2,-1,jvel1);
        const double fy2 = cface(p,a,b,2,0,jvel2);
        dy = (jvel2*fy2 - jvel1*fy1)/(p->DXM);
    }

    const double fz1 = cface(p,a,b,3,-1,kvel1);
    const double fz2 = cface(p,a,b,3,0,kvel2);
    const double dz = (kvel2*fz2 - kvel1*fz1)/(p->DXM);

    return -dx-dy-dz;
}

template<typename GenericField>
double hric_mod::cface(lexer *p, fdm *a, const GenericField& b, int dir, int pos, double uwind)
{
    double cj,cj_,cj_s,cj_ss;
    double cc,cc_,cu,cd;
    double umax;
    double gradx,grady,gradz;

    if(dir==1)
    {
        if(uwind>=0.0)
        {
            cc = b(i+pos,j,k);
            cu = b(i-1+pos,j,k);
            cd = b(i+1+pos,j,k);
            umax = 0.5*(a->u(i+pos,j,k)+a->u(i-1+pos,j,k));
        }
        else
        {
            cc = b(i+1+pos,j,k);
            cu = b(i+2+pos,j,k);
            cd = b(i-0+pos,j,k);
            umax = 0.5*(a->u(i+1-pos,j,k)+a->u(i-0-pos,j,k));
        }
    }
    else if(dir==2)
    {
        if(uwind>=0.0)
        {
            cc = b(i,j+pos,k);
            cu = b(i,j-1+pos,k);
            cd = b(i,j+1+pos,k);
            umax = 0.5*(a->v(i,j+pos,k)+a->v(i,j-1+pos,k));
        }
        else
        {
            cc = b(i,j+1+pos,k);
            cu = b(i,j+2+pos,k);
            cd = b(i,j-0+pos,k);
            umax = 0.5*(a->v(i,j+1-pos,k)+a->v(i,j-0-pos,k));
        }
    }
    else if(dir==3)
    {
        if(uwind>=0.0)
        {
            cc = b(i,j,k+0+pos);
            cu = b(i,j,k-1+pos);
            cd = b(i,j,k+1+pos);
            umax = 0.5*(a->w(i,j,k+pos)+a->w(i,j,k-1+pos));
        }
        else if(uwind<0.0)
        {
            cc = b(i,j,k+1+pos);
            cu = b(i,j,k+2+pos);
            cd = b(i,j,k-0+pos);
            umax = 0.5*(a->w(i,j,k+1-pos)+a->w(i,j,k-0-pos));
        }
    }

    cc_ = (cc-cu)/(fabs(cd-cu)>1.0e-20?(cd-cu):1.0e20);

    if(cc_<0.0)
    cj_ = cc_;
    else if(cc_>=0.0 && cc_<0.5)
    cj_ = 2.0*cc_;
    else if(cc_>=0.5 && cc_<1.0)
    cj_ = 1.0;
    else if(cc_>=1.0)
    cj_ = cc_;

    if(cc_>=0.0 && cc_<1.0)
    cj_s = std::min(cj_,(6.0*cc_+3.0)/8.0);
    else
    cj_s = cc_;

    if(uwind>=0.0)
    {
        gradx = fabs((b(i+1+pos,j,k)-b(i-1+pos,j,k))/(2.0*p->DXM));
        grady = fabs((b(i,j+1+pos,k)-b(i,j-1+pos,k))/(2.0*p->DXM));
        gradz = fabs((b(i,j,k+1+pos)-b(i,j,k-1+pos))/(2.0*p->DXM));
    }
    else
    {
        gradx = fabs((b(i+1+1+pos,j,k)-b(i-1+1+pos,j,k))/(2.0*p->DXM));
        grady = fabs((b(i,j+1+1+pos,k)-b(i,j-1+1+pos,k))/(2.0*p->DXM));
        gradz = fabs((b(i,j,k+1+1+pos)-b(i,j,k-1+1+pos))/(2.0*p->DXM));
    }

    const double vl = sqrt(gradx*gradx + grady*grady + gradz*gradz);

    double costheta = 0.0;
    if(dir==1)
    costheta = gradx/(vl>1.0e-20?vl:1.0e20);
    if(dir==2)
    costheta = grady/(vl>1.0e-20?vl:1.0e20);
    if(dir==3)
    costheta = gradz/(vl>1.0e-20?vl:1.0e20);

    cj_ss = cj_*sqrt(costheta) + cj_s*(1.0-sqrt(costheta));

    return cj_ss*(cd-cu) + cu;
}
