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

#include"weno_hj_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"
#include"flux_HJ_CDS2_2D.h"
#include"flux_HJ_CDS2_vrans_2D.h"

weno_hj_nug::weno_hj_nug(lexer* p) : weno_nug_func(p)
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

void weno_hj_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;

    if(ipol==1)
    {
        uf=1;
        FIELDLOOP_INC_MEMBER(a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXN.data(),p->DYP.data(),p->DZP.data());
        )
    }
    else if(ipol==2 && p->j_dir==1)
    {
        vf=1;
        FIELDLOOP_INC_MEMBER(a,G,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZP.data());
        )
    }
    else if(ipol==3)
    {
        wf=1;
        FIELDLOOP_INC_MEMBER(a,H,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXP.data(),p->DYP.data(),p->DZN.data());
        )
    }
    else if(ipol==4)
    {
        FIELDLOOP_INC_MEMBER(a,L,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXP.data(),p->DYP.data(),p->DZP.data());
        )
    }
}

template<typename GenericField>
inline double weno_hj_nug::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel, double* DX, double* DY, double* DZ)
{
    double iadvec, ivel2, jadvec, jvel2, kadvec, kvel2;

    pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
    pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
    pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

    double L = -iadvec*fx(p,a,b,DX,iadvec);

    if(p->j_dir==1)
    L -= jadvec*fy(p,a,b,DY,jadvec);

    L -= kadvec*fz(p,a,b,DZ,kadvec);

    return L;
}

template<typename GenericField>
double weno_hj_nug::fx(lexer* p, fdm* a, const GenericField& b, const double* DX, double advec)
{
    if(advec==0.0) return 0.0;

    if(advec>0.0)
    {
        q1 = (b(i-2,j,k) - b(i-3,j,k))/DX[IM3];
        q2 = (b(i-1,j,k) - b(i-2,j,k))/DX[IM2];
        q3 = (b(i,j,k)   - b(i-1,j,k))/DX[IM1];
        q4 = (b(i+1,j,k) - b(i,j,k)  )/DX[IP];
        q5 = (b(i+2,j,k) - b(i+1,j,k))/DX[IP1];
        is_min_x();
        weight_min_x();

        return w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
             + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
    }
    else
    {
        q1 = (b(i-1,j,k) - b(i-2,j,k))/DX[IM2];
        q2 = (b(i,j,k)   - b(i-1,j,k))/DX[IM1];
        q3 = (b(i+1,j,k) - b(i,j,k)  )/DX[IP];
        q4 = (b(i+2,j,k) - b(i+1,j,k))/DX[IP1];
        q5 = (b(i+3,j,k) - b(i+2,j,k))/DX[IP2];
        is_max_x();
        weight_max_x();

        return w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
             + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
    }
}

template<typename GenericField>
double weno_hj_nug::fy(lexer* p, fdm* a, const GenericField& b, const double* DY, double advec)
{
    if(advec==0.0) return 0.0;

    if(advec>0.0)
    {
        q1 = (b(i,j-2,k) - b(i,j-3,k))/DY[JM3];
        q2 = (b(i,j-1,k) - b(i,j-2,k))/DY[JM2];
        q3 = (b(i,j,k)   - b(i,j-1,k))/DY[JM1];
        q4 = (b(i,j+1,k) - b(i,j,k)  )/DY[JP];
        q5 = (b(i,j+2,k) - b(i,j+1,k))/DY[JP1];
        is_min_y();
        weight_min_y();

        return w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
             + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
    }
    else
    {
        q1 = (b(i,j-1,k) - b(i,j-2,k))/DY[JM2];
        q2 = (b(i,j,k)   - b(i,j-1,k))/DY[JM1];
        q3 = (b(i,j+1,k) - b(i,j,k)  )/DY[JP];
        q4 = (b(i,j+2,k) - b(i,j+1,k))/DY[JP1];
        q5 = (b(i,j+3,k) - b(i,j+2,k))/DY[JP2];
        is_max_y();
        weight_max_y();

        return w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
             + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
    }
}

template<typename GenericField>
double weno_hj_nug::fz(lexer* p, fdm* a, const GenericField& b, const double* DZ, double advec)
{
    if(advec==0.0) return 0.0;

    if(advec>0.0)
    {
        q1 = (b(i,j,k-2) - b(i,j,k-3))/DZ[KM3];
        q2 = (b(i,j,k-1) - b(i,j,k-2))/DZ[KM2];
        q3 = (b(i,j,k)   - b(i,j,k-1))/DZ[KM1];
        q4 = (b(i,j,k+1) - b(i,j,k)  )/DZ[KP];
        q5 = (b(i,j,k+2) - b(i,j,k+1))/DZ[KP1];
        is_min_z();
        weight_min_z();

        return w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
             + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
    }
    else
    {
        q1 = (b(i,j,k-1) - b(i,j,k-2))/DZ[KM2];
        q2 = (b(i,j,k)   - b(i,j,k-1))/DZ[KM1];
        q3 = (b(i,j,k+1) - b(i,j,k)  )/DZ[KP];
        q4 = (b(i,j,k+2) - b(i,j,k+1))/DZ[KP1];
        q5 = (b(i,j,k+3) - b(i,j,k+2))/DZ[KP2];
        is_max_z();
        weight_max_z();

        return w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
             + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
             + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
    }
}
