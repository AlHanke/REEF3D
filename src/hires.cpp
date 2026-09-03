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

#include"hires.h"
#include"lexer.h"
#include"fdm.h"
#include"minmod.h"
#include"vanleer.h"
#include"umist.h"
#include"vanalbada.h"
#include"superbee.h"
#include"smart.h"
#include"limo3.h"
#include"tvdvof.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"

hires::hires (lexer *p, int limiter)
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


    if(limiter==10)
    plim = new minmod(p);

    else if(limiter==11)
    plim = new vanleer(p);

    else if(limiter==12)
    plim = new umist(p);

    else if(limiter==13)
    plim = new vanalbada(p);

    else if(limiter==14)
    plim = new superbee(p);

    else if(limiter==15)
    plim = new smart(p);

    else if(limiter==16)
    plim = new limo3(p);

    else if(limiter==42)
    plim = new tvdvof(p);
}

void hires::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    {
        FIELDLOOP_INC_MEMBER(a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZN.data());
        )
    }
    else if(ipol==2 && p->j_dir==1)
    {
        FIELDLOOP_INC_MEMBER(a,G,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN.data(),p->DYP.data(),p->DZN.data());
        )
    }
    else if(ipol==3)
    {
        FIELDLOOP_INC_MEMBER(a,H,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZP.data());
        )
    }
    else if(ipol==4)
    {
        FIELDLOOP_INC_MEMBER(a,L,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN.data(),p->DYN.data(),p->DZN.data());
        )
    }
}

template<typename GenericField>
double hires::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel, double *DX,double *DY, double *DZ)
{
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    double udir,vdir,wdir;
    udir=vdir=wdir=0.0;

    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    // x-dir
    if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;

    double dx = udir*(ivel2*(b(i,j,k) + 0.5*plim->iphi(b,0,-1,1,0)*(b(i+1,j,k)-b(i,j,k)))
              - ivel1*(b(i-1,j,k) + 0.5*plim->iphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i-1,j,k))))/DX[IM1]
              + (1.0-udir)*(ivel2*(b(i+1,j,k) - 0.5*plim->iphi(b,1,0,2,1)*(b(i+2,j,k)-b(i+1,j,k)))
              - ivel1*(b(i,j,k) - 0.5*plim->iphi(b,0,-1,1,0)*(b(i+1,j,k)-b(i,j,k))))/DX[IP];

    // y-dir
    double dy=0.0;
    if(p->j_dir==1)
    {
        if(0.5*(jvel1+jvel2)>=0.0)
        vdir=1.0;

        dy = vdir*(jvel2*(b(i,j,k) + 0.5*plim->jphi(b,0,-1,1,0)*(b(i,j+1,k)-b(i,j,k)))
           - jvel1*(b(i,j-1,k) + 0.5*plim->jphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i,j-1,k))))/DY[JM1]
           + (1.0-vdir)*(jvel2*(b(i,j+1,k) - 0.5*plim->jphi(b,1,0,2,1)*(b(i,j+2,k)-b(i,j+1,k)))
           - jvel1*(b(i,j,k) - 0.5*plim->jphi(b,0,-1,1,0)*(b(i,j+1,k)-b(i,j,k))))/DY[JP];
    }

    // z-dir
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;

    double dz = wdir*(kvel2*(b(i,j,k) + 0.5*plim->kphi(b,0,-1,1,0)*(b(i,j,k+1)-b(i,j,k)))
              -  kvel1*(b(i,j,k-1) + 0.5*plim->kphi(b,-1,-2,0,-1)*(b(i,j,k)-b(i,j,k-1))))/DZ[KM1]
              + (1.0-wdir)*(kvel2*(b(i,j,k+1) - 0.5*plim->kphi(b,1,0,2,1)*(b(i,j,k+2)-b(i,j,k+1)))
              - kvel1*(b(i,j,k) - 0.5*plim->kphi(b,0,-1,1,0)*(b(i,j,k+1)-b(i,j,k))))/DZ[KP];

    return -dx-dy-dz;
}
