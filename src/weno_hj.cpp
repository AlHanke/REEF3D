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

#include"weno_hj.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"
#include"flux_HJ_CDS2_2D.h"
#include"flux_HJ_CDS2_vrans_2D.h"

weno_hj::weno_hj(lexer* p)
{
    if(p->j_dir==0)
    {
        if(p->B269>=1 || p->S10==2)
        pflux = new flux_HJ_CDS2_vrans_2D;

        else if(p->B269==0 && p->S10!=2)
        pflux = new flux_HJ_CDS2_2D;
    }
    else if(p->j_dir==1)
    {
        if(p->B269>=1 || p->S10==2)
        pflux = new flux_HJ_CDS2_vrans;

        else if(p->B269==0 && p->S10!=2)
        pflux = new flux_HJ_CDS2;
    }
}

void weno_hj::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    {
        n=0;
        FIELDLOOP_INC_MEMBER(a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DRDXN,p->DSDYP,p->DTDZP); ++n;
        )
    }
    else if(ipol==2 && p->j_dir==1)
    {
        n=0;
        FIELDLOOP_INC_MEMBER(a,G,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DRDXP,p->DSDYN,p->DTDZP); ++n;
        )
    }
    else if(ipol==3)
    {
        n=0;
        FIELDLOOP_INC_MEMBER(a,H,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZN); ++n;
        )
    }
    else if(ipol==4)
    {
        n=0;
        FIELDLOOP_INC_MEMBER(a,L,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DRDXP,p->DSDYP,p->DTDZP); ++n;
        )
    }
}

template<typename GenericField>
inline double weno_hj::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel, double* DRDX, double* DSDY, double* DTDZ)
{
    double iadvec, ivel2, jadvec, jvel2, kadvec, kvel2;

    pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
    pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
    pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

    return - iadvec*ddx(p,a,b,iadvec)*DRDX[IP]
           - jadvec*ddy(p,a,b,jadvec)*DSDY[JP]
           - kadvec*ddz(p,a,b,kadvec)*DTDZ[KP];
}

template<typename GenericField>
inline double weno_hj::ddx(lexer* p, fdm* a, const GenericField& b, double advec)
{
    if(advec==0.0) return 0.0;

    double q1,q2,q3,q4,q5;

    if(advec>0.0)
    {
        q1 = (b(i-2,j,k) - b(i-3,j,k))/p->DRM;
        q2 = (b(i-1,j,k) - b(i-2,j,k))/p->DRM;
        q3 = (b(i,j,k)   - b(i-1,j,k))/p->DRM;
        q4 = (b(i+1,j,k) - b(i,j,k)  )/p->DRM;
        q5 = (b(i+2,j,k) - b(i+1,j,k))/p->DRM;
    }
    else
    {
        q1 = (b(i+3,j,k) - b(i+2,j,k))/p->DRM;
        q2 = (b(i+2,j,k) - b(i+1,j,k))/p->DRM;
        q3 = (b(i+1,j,k) - b(i,j,k)  )/p->DRM;
        q4 = (b(i,j,k)   - b(i-1,j,k))/p->DRM;
        q5 = (b(i-1,j,k) - b(i-2,j,k))/p->DRM;
    }

    const double is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
    const double is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
    const double is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
    const double a1 = tenth/pow(epsilon+is1, 2.0);
    const double a2 = sixten/pow(epsilon+is2, 2.0);
    const double a3 = treten/pow(epsilon+is3, 2.0);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}

template<typename GenericField>
inline double weno_hj::ddy(lexer* p, fdm* a, const GenericField& b, double advec)
{
    if(advec==0.0) return 0.0;

    double q1,q2,q3,q4,q5;

    if(advec>0.0)
    {
        q1 = (b(i,j-2,k) - b(i,j-3,k))/p->DSM;
        q2 = (b(i,j-1,k) - b(i,j-2,k))/p->DSM;
        q3 = (b(i,j,k)   - b(i,j-1,k))/p->DSM;
        q4 = (b(i,j+1,k) - b(i,j,k)  )/p->DSM;
        q5 = (b(i,j+2,k) - b(i,j+1,k))/p->DSM;
    }
    else
    {
        q1 = (b(i,j+3,k) - b(i,j+2,k))/p->DSM;
        q2 = (b(i,j+2,k) - b(i,j+1,k))/p->DSM;
        q3 = (b(i,j+1,k) - b(i,j,k)  )/p->DSM;
        q4 = (b(i,j,k)   - b(i,j-1,k))/p->DSM;
        q5 = (b(i,j-1,k) - b(i,j-2,k))/p->DSM;
    }

    const double is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
    const double is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
    const double is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
    const double a1 = tenth/pow(epsilon+is1, 2.0);
    const double a2 = sixten/pow(epsilon+is2, 2.0);
    const double a3 = treten/pow(epsilon+is3, 2.0);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}

template<typename GenericField>
inline double weno_hj::ddz(lexer* p, fdm* a, const GenericField& b, double advec)
{
    if(advec==0.0) return 0.0;

    double q1,q2,q3,q4,q5;

    if(advec>0.0)
    {
        q1 = (b(i,j,k-2) - b(i,j,k-3))/p->DTM;
        q2 = (b(i,j,k-1) - b(i,j,k-2))/p->DTM;
        q3 = (b(i,j,k)   - b(i,j,k-1))/p->DTM;
        q4 = (b(i,j,k+1) - b(i,j,k)  )/p->DTM;
        q5 = (b(i,j,k+2) - b(i,j,k+1))/p->DTM;
    }
    else
    {
        q1 = (b(i,j,k+3) - b(i,j,k+2))/p->DTM;
        q2 = (b(i,j,k+2) - b(i,j,k+1))/p->DTM;
        q3 = (b(i,j,k+1) - b(i,j,k)  )/p->DTM;
        q4 = (b(i,j,k)   - b(i,j,k-1))/p->DTM;
        q5 = (b(i,j,k-1) - b(i,j,k-2))/p->DTM;
    }

    const double is1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
    const double is2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
    const double is3 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
    const double a1 = tenth/pow(epsilon+is1, 2.0);
    const double a2 = sixten/pow(epsilon+is2, 2.0);
    const double a3 = treten/pow(epsilon+is3, 2.0);
    const double asum = a1+a2+a3;
    const double w1=a1/asum, w2=a2/asum, w3=a3/asum;

    return w1*( q1*third - q2*sevsix + q3*elvsix)
         + w2*(-q2*sixth + q3*fivsix + q4*third)
         + w3*( q3*third + q4*fivsix - q5*sixth);
}
