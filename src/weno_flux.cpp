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

#include"weno_flux.h"
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

weno_flux::weno_flux(lexer* p)
{
    if(p->j_dir==0)
    {
        if(p->B200>=1 || p->S10==2)
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
        if(p->B200>=1 || p->S10==2)
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

inline void weno_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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
inline double weno_flux::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    const double fu1 = fx(p,b,uvel,ipol,ivel1,-1);
    const double fu2 = fx(p,b,uvel,ipol,ivel2,0);

    double fv1=0.0,fv2=0.0;
    if(p->j_dir==1)
    {
        fv1 = fy(p,b,vvel,ipol,jvel1,-1);
        fv2 = fy(p,b,vvel,ipol,jvel2,0);
    }

    const double fw1 = fz(p,b,wvel,ipol,kvel1,-1);
    const double fw2 = fz(p,b,wvel,ipol,kvel2,0);

    return - ((ivel2*fu2-ivel1*fu1)/p->DXM)
            - ((jvel2*fv2-jvel1*fv1)/p->DYM)
            - ((kvel2*fw2-kvel1*fw1)/p->DZM);
}

template<typename GenericField>
inline double weno_flux::fx(lexer *p, const GenericField& b, const GenericField& uvel, int ipol, double advec, int di)
{
    double q1,q2,q3,q4,q5;

    if(advec==0.0) return 0.0;
    else if(advec>0.0)
    {
        q1 = b(i+di-2,j,k);
        q2 = b(i+di-1,j,k);
        q3 = b(i+di,j,k);
        q4 = b(i+di+1,j,k);
        q5 = b(i+di+2,j,k);
    }
    else if(advec<0.0)
    {
        q1 = b(i+di+3,j,k);
        q2 = b(i+di+2,j,k);
        q3 = b(i+di+1,j,k);
        q4 = b(i+di,j,k);
        q5 = b(i+di-1,j,k);
    }

    const double q123a = q1 - 2.0*q2 + q3;
    const double q123b = q1 - 4.0*q2 + 3.0*q3;
    const double q234a = q2 - 2.0*q3 + q4;
    const double q234b = q2 - q4;
    const double q345a = q3 - 2.0*q4 + q5;
    const double q345b = 3.0*q3 - 4.0*q4 + q5;

    const double is1 = tttw*q123a*q123a + fourth*q123b*q123b + epsilon;
    const double is2 = tttw*q234a*q234a + fourth*q234b*q234b + epsilon;
    const double is3 = tttw*q345a*q345a + fourth*q345b*q345b + epsilon;
    const double a1 = tenth/(is1*is1);
    const double a2 = sixten/(is2*is2);
    const double a3 = treten/(is3*is3);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}

template<typename GenericField>
inline double weno_flux::fy(lexer *p, const GenericField& b, const GenericField& vvel, int ipol, double advec, int dj)
{
    double q1,q2,q3,q4,q5;

    if(advec==0.0) return 0.0;
    else if(advec>0.0)
    {
        q1 = b(i,j+dj-2,k);
        q2 = b(i,j+dj-1,k);
        q3 = b(i,j+dj,k);
        q4 = b(i,j+dj+1,k);
        q5 = b(i,j+dj+2,k);
    }
    else if(advec<0.0)
    {
        q1 = b(i,j+dj+3,k);
        q2 = b(i,j+dj+2,k);
        q3 = b(i,j+dj+1,k);
        q4 = b(i,j+dj,k);
        q5 = b(i,j+dj-1,k);
    }

    const double q123a = q1 - 2.0*q2 + q3;
    const double q123b = q1 - 4.0*q2 + 3.0*q3;
    const double q234a = q2 - 2.0*q3 + q4;
    const double q234b = q2 - q4;
    const double q345a = q3 - 2.0*q4 + q5;
    const double q345b = 3.0*q3 - 4.0*q4 + q5;

    const double is1 = tttw*q123a*q123a + fourth*q123b*q123b + epsilon;
    const double is2 = tttw*q234a*q234a + fourth*q234b*q234b + epsilon;
    const double is3 = tttw*q345a*q345a + fourth*q345b*q345b + epsilon;
    const double a1 = tenth/(is1*is1);
    const double a2 = sixten/(is2*is2);
    const double a3 = treten/(is3*is3);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}

template<typename GenericField>
inline double weno_flux::fz(lexer *p, const GenericField& b, const GenericField& wvel, int ipol, double advec, int dk)
{
    double q1,q2,q3,q4,q5;

    if(advec==0.0) return 0.0;
    else if(advec>0.0)
    {
        q1 = b(i,j,k+dk-2);
        q2 = b(i,j,k+dk-1);
        q3 = b(i,j,k+dk);
        q4 = b(i,j,k+dk+1);
        q5 = b(i,j,k+dk+2);
    }
    else if(advec<0.0)
    {
        q1 = b(i,j,k+dk+3);
        q2 = b(i,j,k+dk+2);
        q3 = b(i,j,k+dk+1);
        q4 = b(i,j,k+dk);
        q5 = b(i,j,k+dk-1);
    }

    const double q123a = q1 - 2.0*q2 + q3;
    const double q123b = q1 - 4.0*q2 + 3.0*q3;
    const double q234a = q2 - 2.0*q3 + q4;
    const double q234b = q2 - q4;
    const double q345a = q3 - 2.0*q4 + q5;
    const double q345b = 3.0*q3 - 4.0*q4 + q5;

    const double is1 = tttw*q123a*q123a + fourth*q123b*q123b + epsilon;
    const double is2 = tttw*q234a*q234a + fourth*q234b*q234b + epsilon;
    const double is3 = tttw*q345a*q345a + fourth*q345b*q345b + epsilon;
    const double a1 = tenth/(is1*is1);
    const double a2 = sixten/(is2*is2);
    const double a3 = treten/(is3*is3);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}
