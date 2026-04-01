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

#include"reinidisc_f.h"
#include"lexer.h"
#include"field.h"

void reinidisc_f::start(lexer *p, fdm*, ghostcell *pgc, field &f, field &L, int ipol) noexcept
{
    if(ipol==4 || ipol==5)
    {
        const bool is3D = (p->j_dir == 1);

        L.setVal(0.0);

        FIELDLOOP_INC(
            L,
            FIELD_CONST_INC(f),
            L(i,j,k) = disc(p,f,is3D);
        )
    }
}

template<typename GenericFieldConst>
inline double reinidisc_f::disc(lexer *p, const GenericFieldConst &f, const bool is3D) noexcept
{
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    const double lsv = f(i,j,k);
    const double lsv2 = lsv*lsv;
    const double lsSig = (fabs(lsv) < 1.0e-8) ? 1.0 : std::copysign(1.0, lsv);

// x
    const double xmin = (lsv-f(i-1,j,k))/p->DXP[IM1];
    const double xplus = (f(i+1,j,k)-lsv)/p->DXP[IP];
    const double xmin_s = xmin * lsSig;
    const double xplus_s = xplus * lsSig;

    if     (xmin_s>0.0 && xplus_s>-xmin_s) dx = ddwenox(f,1.0);
    else if(xplus_s<0.0 && xmin_s<-xplus_s) dx = ddwenox(f,-1.0);

// y
    if(is3D)
    {
        const double ymin = (lsv-f(i,j-1,k))/p->DYP[JM1];
        const double yplus = (f(i,j+1,k)-lsv)/p->DYP[JP];
        const double ymin_s = ymin * lsSig;
        const double yplus_s = yplus * lsSig;

        if     (ymin_s>0.0 && yplus_s>-ymin_s) dy = ddwenoy(f,1.0);
        else if(yplus_s<0.0 && ymin_s<-yplus_s) dy = ddwenoy(f,-1.0);
    }

// z
    const double zmin = (lsv-f(i,j,k-1))/p->DZP[KM1];
    const double zplus = (f(i,j,k+1)-lsv)/p->DZP[KP];
    const double zmin_s = zmin * lsSig;
    const double zplus_s = zplus * lsSig;

    if     (zmin_s>0.0 && zplus_s>-zmin_s) dz = ddwenoz(f,1.0);
    else if(zplus_s<0.0 && zmin_s<-zplus_s) dz = ddwenoz(f,-1.0);

    const double dnorm = sqrt(dx*dx + dy*dy + dz*dz);

    const double deltax = (is3D) ? (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]) : (1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);

    double sign = lsv/sqrt(lsv2 + dnorm*dnorm*deltax*deltax);
    if(std::isnan(sign))
    sign = 1.0;

    return sign * (1.0 - dnorm);
}
