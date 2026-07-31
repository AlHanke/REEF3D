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

#include"cds2_alt.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"

cds2_alt::cds2_alt(lexer *p)
{
    if(p->B200>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans;

    else
    pflux = new flux_HJ_CDS2;
}

void cds2_alt::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    {
        FIELDLOOP_INC_MEMBER(a,F,
            FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_INC(vvel); FIELD_CONST_INC(wvel),
            F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP.data(),p->DYN.data(),p->DZN.data());
        )
    }
    else if(ipol==2)
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
double cds2_alt::aij(lexer* p, fdm* a, const GenericField& b, int ipol, const GenericField& uvel, const GenericField& vvel, const GenericField& wvel, double *DX, double *DY, double *DZ)
{
    double iadvec, jadvec, kadvec, temp;
    pflux->u_flux(a,ipol,uvel,iadvec,temp);
    pflux->v_flux(a,ipol,vvel,jadvec,temp);
    pflux->w_flux(a,ipol,wvel,kadvec,temp);

    const double dx = iadvec*(0.5*(b(i,j,k) + b(i+1,j,k))  -  0.5*(b(i-1,j,k) +  b(i,j,k)))/DX[IP];

    const double dy = jadvec*(0.5*(b(i,j,k) + b(i,j+1,k))  -  0.5*(b(i,j-1,k) +  b(i,j,k)))/DY[JP];

    const double dz = kadvec*(0.5*(b(i,j,k) + b(i,j,k+1))  -  0.5*(b(i,j,k-1) +  b(i,j,k)))/DZ[KP];

    return -dx-dy-dz;
}
