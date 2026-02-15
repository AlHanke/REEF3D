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
        if(p->B269==0)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_2D(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_2D(p);
        }
        else if(p->B269>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans_2D(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans_2D(p);
        }
    }
    else if(p->j_dir==1)
    {
        if(p->B269==0)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2(p);
        }
        else if(p->B269>=1 || p->S10==2)
        {
            if(p->D11==1)
            pflux = new flux_face_FOU_vrans(p);
            else if(p->D11==2)
            pflux = new flux_face_CDS2_vrans(p);
        }
    }
}

void weno3_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    uf=vf=wf=0;

    if(ipol==1)
    {
        uf=1;
        ULOOP
        a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZN.data());
    }
    else if(ipol==2 && p->j_dir==1)
    {
        vf=1;
        VLOOP
        a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN.data(),p->DYP.data(),p->DZN.data());
    }
    else if(ipol==3)
    {
        wf=1;
        WLOOP
        a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZP.data());
    }
    else if(ipol==4)
    {
        LOOP
        a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZN.data());
    }
}

inline double weno3_flux::aij(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel, double *DX, double *DY, double *DZ)
{
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    i -= 1;
    fu1 = fx(p,b,ivel1);
    i += 1;

    fu2 = fx(p,b,ivel2);

    if(p->j_dir==1)
    {
        j -=1 ;
        fv1 = fy(p,b,jvel1);
        j += 1;

        fv2 = fy(p,b,jvel2);
    }
    else
    {
        fv1 = fv2 = 0.0;
    }

    k -= 1;
    fw1 = fz(p,b,kvel1);
    k += 1;

    fw2 = fz(p,b,kvel2);

    L = - ((ivel2*fu2-ivel1*fu1)/DX[IP])
        - ((jvel2*fv2-jvel1*fv1)/DY[JP])
        - ((kvel2*fw2-kvel1*fw1)/DZ[KP]);

    return L;
}

inline double weno3_flux::fx(lexer *p, field& b, double advec)
{
    grad = 0.0;

    if(advec>0.0)
    {
        q1 = b(i-1,j,k);
        q2 = b(i,j,k);
        q3 = b(i+1,j,k);

        is_min_x();
        weight_min_x();

        grad = w1x*(qfx[IP][uf][0][0]*q2 + qfx[IP][uf][0][1]*q3)
             + w2x*(qfx[IP][uf][1][0]*q2 - qfx[IP][uf][1][1]*q1);
    }
    else if(advec<0.0)
    {
        q1 = b(i,j,k);
        q2 = b(i+1,j,k);
        q3 = b(i+2,j,k);

        is_max_x();
        weight_max_x();

        grad = w1x*(qfx[IP][uf][2][0]*q2 - qfx[IP][uf][2][1]*q3)
             + w2x*(qfx[IP][uf][3][0]*q1 + qfx[IP][uf][3][1]*q2);
    }

    return grad;
}

inline double weno3_flux::fy(lexer *p, field& b, double advec)
{
    grad = 0.0;

    if(advec>0.0)
    {
        q1 = b(i,j-1,k);
        q2 = b(i,j,k);
        q3 = b(i,j+1,k);

        is_min_y();
        weight_min_y();

        grad = w1y*(qfy[JP][vf][0][0]*q2 + qfy[JP][vf][0][1]*q3)
             + w2y*(qfy[JP][vf][1][0]*q2 - qfy[JP][vf][1][1]*q1);
    }
    else if(advec<0.0)
    {
        q1 = b(i,j,k);
        q2 = b(i,j+1,k);
        q3 = b(i,j+2,k);

        is_max_y();
        weight_max_y();

        grad = w1y*(qfy[JP][vf][2][0]*q2 - qfy[JP][vf][2][1]*q3)
             + w2y*(qfy[JP][vf][3][0]*q1 + qfy[JP][vf][3][1]*q2);
    }

    return grad;
}

inline double weno3_flux::fz(lexer *p, field& b, double advec)
{
    grad = 0.0;

    if(advec>0.0)
    {
        q1 = b(i,j,k-1);
        q2 = b(i,j,k);
        q3 = b(i,j,k+1);

        is_min_z();
        weight_min_z();

        grad = w1z*(qfz[KP][wf][0][0]*q2 + qfz[KP][wf][0][1]*q3)
             + w2z*(qfz[KP][wf][1][0]*q2 - qfz[KP][wf][1][1]*q1);
    }
    else if(advec<0.0)
    {
        q1 = b(i,j,k);
        q2 = b(i,j,k+1);
        q3 = b(i,j,k+2);

        is_max_z();
        weight_max_z();

        grad = w1z*(qfz[KP][wf][2][0]*q2 - qfz[KP][wf][2][1]*q3)
             + w2z*(qfz[KP][wf][3][0]*q1 + qfz[KP][wf][3][1]*q2);
    }

    return grad;
}
