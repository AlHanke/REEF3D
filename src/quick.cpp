/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"quick.h"
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

quick::quick(lexer *p)
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

inline void quick::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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
inline double quick::aij(lexer* p,fdm* a,const GenericField& b,int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel)
{
    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;

    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

    const double ul = ivel1>=0.0 ? 1.0 : 0.0;
    const double ur = ivel2>=0.0 ? 1.0 : 0.0;

    const double dx = (ivel2*(ur*(3.0*b(i+1,j,k) + 6.0*b(i,j,k) - b(i-1,j,k))
                        +(1.0-ur)*(3.0*b(i,j,k) + 6.0*b(i+1,j,k) - b(i+2,j,k)))
                    -  ivel1*(ul*(3.0*b(i,j,k) + 6.0*b(i-1,j,k) - b(i-2,j,k))
                        +(1.0-ul)*(3.0*b(i-1,j,k) + 6.0*b(i,j,k) - b(i+1,j,k))))/(8.0*p->DXM);

    double dy = 0.0;
    if(p->j_dir==1)
    {
        const double vl = jvel1>=0.0 ? 1.0 : 0.0;
        const double vr = jvel2>=0.0 ? 1.0 : 0.0;

        dy = (jvel2*(vr*(3.0*b(i,j+1,k) + 6.0*b(i,j,k) - b(i,j-1,k))
                    +(1.0-vr)*(3.0*b(i,j,k) + 6.0*b(i,j+1,k) - b(i,j+2,k)))
            -  jvel1*(vl*(3.0*b(i,j,k) + 6.0*b(i,j-1,k) - b(i,j-2,k))
                    +(1.0-vl)*(3.0*b(i,j-1,k) + 6.0*b(i,j,k) - b(i,j+1,k))))/(8.0*p->DXM);
    }

    const double wl = kvel1>=0.0 ? 1.0 : 0.0;
    const double wr = kvel2>=0.0 ? 1.0 : 0.0;

    const double dz = (kvel2*(wr*(3.0*b(i,j,k+1) + 6.0*b(i,j,k) - b(i,j,k-1))
                        +(1.0-wr)*(3.0*b(i,j,k) + 6.0*b(i,j,k+1) - b(i,j,k+2)))
                    -  kvel1*(wl*(3.0*b(i,j,k) + 6.0*b(i,j,k-1) - b(i,j,k-2))
                        +(1.0-wl)*(3.0*b(i,j,k-1) + 6.0*b(i,j,k) - b(i,j,k+1))))/(8.0*p->DXM);

    return -dx-dy-dz;
}
