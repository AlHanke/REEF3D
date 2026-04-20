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

#include"iweno_hj_df_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"
#include"flux_HJ_CDS2_2D.h"
#include"flux_HJ_CDS2_vrans_2D.h"

iweno_hj_df_nug::iweno_hj_df_nug(lexer *p)
            :weno_nug_func(p), tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
            sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
            sixten(6.0/10.0),treten(3.0/10.0),epsi(1.0e-6),deltin (1.0/p->DXM)
{
    if(p->j_dir==0)
    {
        if(p->B269>=1 || p->S10==2)
        pflux = new flux_HJ_CDS2_vrans_2D;
        else
        pflux = new flux_HJ_CDS2_2D;
    }
    else if(p->j_dir==1)
    {
        if(p->B269>=1 || p->S10==2)
        pflux = new flux_HJ_CDS2_vrans;
        else
        pflux = new flux_HJ_CDS2;
    }
}

void iweno_hj_df_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;

    if(ipol==1)
    wenoloop1(p,a,b,ipol,uvel,vvel,wvel);

    else if(ipol==2 && p->j_dir==1)
    wenoloop2(p,a,b,ipol,uvel,vvel,wvel);

    else if(ipol==3)
    wenoloop3(p,a,b,ipol,uvel,vvel,wvel);

    else if(ipol==4)
    wenoloop4(p,a,b,ipol,uvel,vvel,wvel);
}

template<typename GenericField>
void iweno_hj_df_nug::wenoloop1(lexer *p, fdm *a, const GenericField& f, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    DX=p->DXN.data();
    DY=p->DYP.data();
    DZ=p->DZP.data();

    uf=1;
    count=0;

    FIELDLOOP_INC_MEMBER(a, F,
        FIELD_CONST_INC(f); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
        {
            pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
            pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
            pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

            if(iadvec>=0.0) { iqmin(p,a,f); is_min_x(); weight_min_x(); aij_south(p,a,f,F); }
            if(iadvec<0.0)  { iqmax(p,a,f); is_max_x(); weight_max_x(); aij_north(p,a,f,F); }

            if(jadvec>=0.0 && p->j_dir==1) { jqmin(p,a,f); is_min_y(); weight_min_y(); aij_east(p,a,f,F); }
            if(jadvec<0.0 && p->j_dir==1)  { jqmax(p,a,f); is_max_y(); weight_max_y(); aij_west(p,a,f,F); }

            if(kadvec>=0.0) { kqmin(p,a,f); is_min_z(); weight_min_z(); aij_bottom(p,a,f,F); }
            if(kadvec<0.0)  { kqmax(p,a,f); is_max_z(); weight_max_z(); aij_top(p,a,f,F); }
            ++count;
        }
    )
}

template<typename GenericField>
void iweno_hj_df_nug::wenoloop2(lexer *p, fdm *a, const GenericField& f, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    DX=p->DXP.data();
    DY=p->DYN.data();
    DZ=p->DZP.data();

    vf=1;
    count=0;

    FIELDLOOP_INC_MEMBER(a, G,
        FIELD_CONST_INC(f); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
        {
            pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
            pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
            pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

            if(iadvec>=0.0) { iqmin(p,a,f); is_min_x(); weight_min_x(); aij_south(p,a,f,G); }
            if(iadvec<0.0)  { iqmax(p,a,f); is_max_x(); weight_max_x(); aij_north(p,a,f,G); }

            if(jadvec>=0.0 && p->j_dir==1) { jqmin(p,a,f); is_min_y(); weight_min_y(); aij_east(p,a,f,G); }
            if(jadvec<0.0 && p->j_dir==1)  { jqmax(p,a,f); is_max_y(); weight_max_y(); aij_west(p,a,f,G); }

            if(kadvec>=0.0) { kqmin(p,a,f); is_min_z(); weight_min_z(); aij_bottom(p,a,f,G); }
            if(kadvec<0.0)  { kqmax(p,a,f); is_max_z(); weight_max_z(); aij_top(p,a,f,G); }
            ++count;
        }
    )
}

template<typename GenericField>
void iweno_hj_df_nug::wenoloop3(lexer *p, fdm *a, const GenericField& f, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    DX=p->DXP.data();
    DY=p->DYP.data();
    DZ=p->DZN.data();

    wf=1;
    count=0;

    FIELDLOOP_INC_MEMBER(a, H,
        FIELD_CONST_INC(f); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
        {
            pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
            pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
            pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

            if(iadvec>=0.0) { iqmin(p,a,f); is_min_x(); weight_min_x(); aij_south(p,a,f,H); }
            if(iadvec<0.0)  { iqmax(p,a,f); is_max_x(); weight_max_x(); aij_north(p,a,f,H); }

            if(jadvec>=0.0 && p->j_dir==1) { jqmin(p,a,f); is_min_y(); weight_min_y(); aij_east(p,a,f,H); }
            if(jadvec<0.0 && p->j_dir==1)  { jqmax(p,a,f); is_max_y(); weight_max_y(); aij_west(p,a,f,H); }

            if(kadvec>=0.0) { kqmin(p,a,f); is_min_z(); weight_min_z(); aij_bottom(p,a,f,H); }
            if(kadvec<0.0)  { kqmax(p,a,f); is_max_z(); weight_max_z(); aij_top(p,a,f,H); }
            ++count;
        }
    )
}

template<typename GenericField>
void iweno_hj_df_nug::wenoloop4(lexer *p, fdm *a, const GenericField& f, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    DX=p->DXP.data();
    DY=p->DYP.data();
    DZ=p->DZP.data();

    uf=vf=wf=0;
    count=0;

    FIELDLOOP_INC_MEMBER(a, L,
        FIELD_CONST_INC(f); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
        {
            pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
            pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
            pflux->w_flux(a,ipol,wvel,kadvec,kvel2);

            if(iadvec>=0.0) { iqmin(p,a,f); is_min_x(); weight_min_x(); aij_south(p,a,f,L); }
            if(iadvec<0.0)  { iqmax(p,a,f); is_max_x(); weight_max_x(); aij_north(p,a,f,L); }

            if(jadvec>=0.0 && p->j_dir==1) { jqmin(p,a,f); is_min_y(); weight_min_y(); aij_east(p,a,f,L); }
            if(jadvec<0.0 && p->j_dir==1)  { jqmax(p,a,f); is_max_y(); weight_max_y(); aij_west(p,a,f,L); }

            if(kadvec>=0.0) { kqmin(p,a,f); is_min_z(); weight_min_z(); aij_bottom(p,a,f,L); }
            if(kadvec<0.0)  { kqmax(p,a,f); is_max_z(); weight_max_z(); aij_top(p,a,f,L); }
            ++count;
        }
    )
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_south(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    -= f(i-3,j,k)*(-w3x*qfx[IP][uf][2][0]/DX[IM3])*iadvec

                + f(i-2,j,k)*(w2x*qfx[IP][uf][1][1]/DX[IM2]
                           - w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2]
                           + w3x*(qfx[IP][uf][2][0])/DX[IM3])*iadvec

                + f(i+2,j,k)*(-w1x*qfx[IP][uf][0][1]/DX[IP1])*iadvec;


    a->M.p[count] = (- w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP]
                    - w2x*(qfx[IP][uf][1][0])/DX[IP]
                    + w1x*(qfx[IP][uf][0][0])/DX[IM1]
                    + w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    + w3x*(qfx[IP][uf][2][1])/DX[IM1])*iadvec;

    a->M.s[count] = (- w1x*qfx[IP][uf][0][0]/DX[IM1]
                    - w2x*(1.0 - qfx[IP][uf][1][0] + qfx[IP][uf][1][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][2][1])/DX[IM1]
                    - w2x*(qfx[IP][uf][1][1])/DX[IM2]
                    + w3x*(1.0 - qfx[IP][uf][2][0] - qfx[IP][uf][2][1])/DX[IM2])*iadvec;

    a->M.n[count] = (  w1x*qfx[IP][uf][0][1]/DX[IP1]
                    + w1x*(1.0 - qfx[IP][uf][0][0] + qfx[IP][uf][0][1])/DX[IP]
                    + w2x*(qfx[IP][uf][1][0])/DX[IP])*iadvec;
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_north(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    -= f(i-2,j,k)*(w3x*qfx[IP][uf][5][1]/DX[IM2])*iadvec

                + f(i+2,j,k)*(-w1x*qfx[IP][uf][3][1]/DX[IP2]
                            + w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                            - w2x*(qfx[IP][uf][4][1])/DX[IP1])*iadvec

                + f(i+3,j,k)*(w1x*qfx[IP][uf][3][1]/DX[IP2])*iadvec;


    a->M.p[count] = (- w1x*(qfx[IP][uf][3][0])/DX[IP]
                    - w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP]
                    - w3x*(qfx[IP][uf][5][0])/DX[IP]
                    + w2x*(qfx[IP][uf][4][0])/DX[IM1]
                    + w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1])*iadvec;

    a->M.s[count] = (- w2x*qfx[IP][uf][4][0]/DX[IM1]
                    - w3x*(1.0 - qfx[IP][uf][5][0] + qfx[IP][uf][5][1])/DX[IM1]
                    - w3x*(qfx[IP][uf][5][1])/DX[IM2])*iadvec;

    a->M.n[count] = (-w1x*(1.0 - qfx[IP][uf][3][0] - qfx[IP][uf][3][1])/DX[IP1]
                    + w2x*(qfx[IP][uf][4][1])/DX[IP1]
                    + w1x*(qfx[IP][uf][3][0])/DX[IP]
                    + w2x*(1.0 - qfx[IP][uf][4][0] + qfx[IP][uf][4][1])/DX[IP]
                    + w3x*qfx[IP][uf][5][0]/DX[IP])*iadvec;
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_east(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    += f(i,j-3,k)*(w3y*qfy[JP][vf][2][0]/DY[JM3])*jadvec

                + f(i,j-2,k)*(-w2y*qfy[JP][vf][1][1]/DY[JM2]
                           + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2]
                           - w3y*(qfy[JP][vf][2][0])/DY[JM3])*jadvec

                + f(i,j+2,k)*(w1y*qfy[JP][vf][0][1]/DY[JP1])*jadvec;


    a->M.p[count] +=(-w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP]
                    - w2y*(qfy[JP][vf][1][0])/DY[JP]
                    + w1y*(qfy[JP][vf][0][0])/DY[JM1]
                    + w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    + w3y*(qfy[JP][vf][2][1])/DY[JM1])*jadvec;

    a->M.e[count] = (-w1y*qfy[JP][vf][0][0]/DY[JM1]
                    - w2y*(1.0 - qfy[JP][vf][1][0] + qfy[JP][vf][1][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][2][1])/DY[JM1]
                    + w2y*(-qfy[JP][vf][1][1])/DY[JM2]
                    + w3y*(1.0 - qfy[JP][vf][2][0] - qfy[JP][vf][2][1])/DY[JM2])*jadvec;

    a->M.w[count] = ( w1y*qfy[JP][vf][0][1]/DY[JP1]
                    + w1y*(1.0 - qfy[JP][vf][0][0] + qfy[JP][vf][0][1])/DY[JP]
                    + w2y*(qfy[JP][vf][1][0])/DY[JP])*jadvec;
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_west(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    += f(i,j-2,k)*(-w3y*qfy[JP][vf][5][1]/DY[JM2])*jadvec

                + f(i,j+2,k)*(w1y*qfy[JP][vf][3][1]/DY[JP2]
                           - w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                           + w2y*(qfy[JP][vf][4][1])/DY[JP1])*jadvec

                + f(i,j+3,k)*(-w1y*qfy[JP][vf][3][1]/DY[JP2])*jadvec;


    a->M.p[count] +=(-w1y*(qfy[JP][vf][3][0])/DY[JP]
                    - w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP]
                    - w3y*(qfy[JP][vf][5][0])/DY[JP]
                    + w2y*(qfy[JP][vf][4][0])/DY[JM1]
                    + w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1])*jadvec;

    a->M.e[count] = (-w2y*qfy[JP][vf][4][0]/DY[JM1]
                    - w3y*(1.0 - qfy[JP][vf][5][0] + qfy[JP][vf][5][1])/DY[JM1]
                    - w3y*(qfy[JP][vf][5][1])/DY[JM2])*jadvec;

    a->M.w[count] = (-w1y*(1.0 - qfy[JP][vf][3][0] - qfy[JP][vf][3][1])/DY[JP1]
                    + w2y*(qfy[JP][vf][4][1])/DY[JP1]
                    + w1y*(qfy[JP][vf][3][0])/DY[JP]
                    + w2y*(1.0 - qfy[JP][vf][4][0] + qfy[JP][vf][4][1])/DY[JP]
                    + w3y*qfy[JP][vf][5][0]/DY[JP])*jadvec;
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_bottom(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    -= f(i,j,k-3)*(-w3z*qfz[KP][wf][2][0]/DZ[KM3])*kadvec

                + f(i,j,k-2)*( w2z*qfz[KP][wf][1][1]/DZ[KM2]
                             - w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2]
                             + w3z*(qfz[KP][wf][2][0])/DZ[KM3])*kadvec

                + f(i,j,k+2)*(-w1z*qfz[KP][wf][0][1]/DZ[KP1])*kadvec;


    a->M.p[count]+= (-w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP]
                    - w2z*(qfz[KP][wf][1][0])/DZ[KP]
                    + w1z*(qfz[KP][wf][0][0])/DZ[KM1]
                    + w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    + w3z*(qfz[KP][wf][2][1])/DZ[KM1])*kadvec;

    a->M.b[count] = (-w1z*qfz[KP][wf][0][0]/DZ[KM1]
                    - w2z*(1.0 - qfz[KP][wf][1][0] + qfz[KP][wf][1][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][2][1])/DZ[KM1]
                    - w2z*(qfz[KP][wf][1][1])/DZ[KM2]
                    + w3z*(1.0 - qfz[KP][wf][2][0] - qfz[KP][wf][2][1])/DZ[KM2])*kadvec;

    a->M.t[count] = ( w1z*qfz[KP][wf][0][1]/DZ[KP1]
                    + w1z*(1.0 - qfz[KP][wf][0][0] + qfz[KP][wf][0][1])/DZ[KP]
                    + w2z*(qfz[KP][wf][1][0])/DZ[KP])*kadvec;
}

template<typename CF, typename MF>
void iweno_hj_df_nug::aij_top(lexer* p, fdm* a, const CF &f, MF &F)
{
    F(i,j,k)    -= f(i,j,k-2)*(w3z*qfz[KP][wf][5][1]/DZ[KM2])*kadvec

                + f(i,j,k+2)*(- w1z*qfz[KP][wf][3][1]/DZ[KP2]
                             + w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                             - w2z*(qfz[KP][wf][4][1])/DZ[KP1])*kadvec

                + f(i,j,k+3)*(w1z*qfz[KP][wf][3][1]/DZ[KP2])*kadvec;


    a->M.p[count]+= (-w1z*(qfz[KP][wf][3][0])/DZ[KP]
                    - w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP]
                    - w3z*(qfz[KP][wf][5][0])/DZ[KP]
                    + w2z*(qfz[KP][wf][4][0])/DZ[KM1]
                    + w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1])*kadvec;

    a->M.b[count] = (-w2z*qfz[KP][wf][4][0]/DZ[KM1]
                    - w3z*(1.0 - qfz[KP][wf][5][0] + qfz[KP][wf][5][1])/DZ[KM1]
                    - w3z*(qfz[KP][wf][5][1])/DZ[KM2])*kadvec;

    a->M.t[count] = (-w1z*(1.0 - qfz[KP][wf][3][0] - qfz[KP][wf][3][1])/DZ[KP1]
                    + w2z*(qfz[KP][wf][4][1])/DZ[KP1]
                    + w1z*(qfz[KP][wf][3][0])/DZ[KP]
                    + w2z*(1.0 - qfz[KP][wf][4][0] + qfz[KP][wf][4][1])/DZ[KP]
                    + w3z*qfz[KP][wf][5][0]/DZ[KP])*kadvec;
}

template<typename GenericField>
void iweno_hj_df_nug::iqmin(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[Im2JK]>0 && p->flagsf4[Im3JK]>0)
    q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];

    if(p->flagsf4[Im1JK]>0 && p->flagsf4[Im2JK]>0)
    q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]>0)
    q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(p->flagsf4[Ip1JK]>0 && p->flagsf4[IJK]>0)
    q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(p->flagsf4[Ip2JK]>0 && p->flagsf4[Ip1JK]>0)
    q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

template<typename GenericField>
void iweno_hj_df_nug::jqmin(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[IJm2K]>0 && p->flagsf4[IJm3K]>0)
    q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];

    if(p->flagsf4[IJm1K]>0 && p->flagsf4[IJm2K]>0)
    q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]>0)
    q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(p->flagsf4[IJp1K]>0 && p->flagsf4[IJK]>0)
    q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(p->flagsf4[IJp2K]>0 && p->flagsf4[IJp1K]>0)
    q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

template<typename GenericField>
void iweno_hj_df_nug::kqmin(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[IJKm2]>0 && p->flagsf4[IJKm3]>0)
    q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];

    if(p->flagsf4[IJKm1]>0 && p->flagsf4[IJKm2]>0)
    q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]>0)
    q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(p->flagsf4[IJKp1]>0 && p->flagsf4[IJKp2]>0)
    q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(p->flagsf4[IJKp2]>0 && p->flagsf4[IJKp1]>0)
    q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

template<typename GenericField>
void iweno_hj_df_nug::iqmax(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[Im1JK]>0 && p->flagsf4[Im1JK]>0)
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[Im1JK]>0)
    q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];

    if(p->flagsf4[Ip1JK]>0 && p->flagsf4[IJK]>0)
    q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];

    if(p->flagsf4[Ip2JK]>0 && p->flagsf4[Ip1JK]>0)
    q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];

    if(p->flagsf4[Ip3JK]>0 && p->flagsf4[Ip2JK]>0)
    q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

template<typename GenericField>
void iweno_hj_df_nug::jqmax(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[IJm1K]>0 && p->flagsf4[IJm2K]>0)
    q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[IJm1K]>0)
    q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];

    if(p->flagsf4[IJp1K]>0 && p->flagsf4[IJK]>0)
    q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];

    if(p->flagsf4[IJp2K]>0 && p->flagsf4[IJp1K]>0)
    q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];

    if(p->flagsf4[IJp3K]>0 && p->flagsf4[IJp2K]>0)
    q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

template<typename GenericField>
void iweno_hj_df_nug::kqmax(lexer *p, fdm *a, const GenericField& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->flagsf4[IJKm1]>0 && p->flagsf4[IJKm2]>0)
    q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];

    if(p->flagsf4[IJK]>0 && p->flagsf4[IJKm1]>0)
    q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];

    if(p->flagsf4[IJKp1]>0 && p->flagsf4[IJK]>0)
    q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];

    if(p->flagsf4[IJKp2]>0 && p->flagsf4[IJKp1]>0)
    q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];

    if(p->flagsf4[IJKp3]>0 && p->flagsf4[IJKp2]>0)
    q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}
