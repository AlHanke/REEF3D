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

#include "reinidisc_f.h"
#include "lexer.h"
#include "field.h"
#include <cstdlib>

reinidisc_f::reinidisc_f(lexer *p) : ddweno_nug_sf(p)
{
    const char* fb = std::getenv("REEF_REINI_FREEZE_BAND");
    freeze_band = (fb != nullptr);
    freeze_fac  = (fb && std::atof(fb) > 0.0) ? std::atof(fb) : 1.0;
}

void reinidisc_f::start(lexer *p, fdm*, ghostcell*, field &f, field &L, int ipol) noexcept
{
    if(ipol==4 || ipol==5)
    {
        L.setVal(0.0);

        if(p->j_dir == 1)
        {
            FIELDLOOP_INC(L, FIELD_CONST_INC(f),
                {
                    #if USE_AMREX
                    if(!_covered_array(i,j,k))
                    #endif
                    L(i,j,k) = disc<true>(p,f);
                }
            )
        }
        else
        {
            FIELDLOOP_INC(L, FIELD_CONST_INC(f),
                {
                    #if USE_AMREX
                    if(!_covered_array(i,j,k))
                    #endif
                    L(i,j,k) = disc<false>(p,f);
                }
            )
        }
    }
}

template<bool Is3D, typename GenericFieldConst>
inline double reinidisc_f::disc(lexer *p, const GenericFieldConst &f) noexcept
{
    // REEF_REINI_FREEZE_BAND: freeze the density band during per-step reinit. density =
    // heaviside_ls(phi,psi) depends on phi only for |phi|<psi (heaviside_ls saturates outside),
    // so reinit restoring |grad phi|=1 inside the band shifts rho every step and closes a
    // velocity<->reinit feedback loop (tighter reinit made umax WORSE, not better). Returning 0
    // leaves band cells un-reinitialised, so rho follows the (barely-moving) advected phi and the
    // loop opens. count>0 only: the count==0 SDF initialisation still needs the full band reinit
    // (F40=0 leaves init broken). Caveat: band phi is no longer a clean SDF, so any curvature/
    // surface-tension consumer of band phi loses accuracy -- acceptable for the at-rest/mild case
    // under test; the robust variant is a |grad phi|-normalised Heaviside in density_f::roface.
    if(freeze_band && p->count>0 && std::fabs(f(i,j,k)) < freeze_fac*p->psi)
    return 0.0;

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    const double lsv = f(i,j,k);
    const double lsv2 = lsv*lsv;
    const double lsSig = (lsv2 < 1.0e-16) ? 1.0 : std::copysign(1.0, lsv);

// x
    const double xmin = (lsv-f(i-1,j,k))/p->DXP[IM1];
    const double xplus = (f(i+1,j,k)-lsv)/p->DXP[IP];
    const double xmin_s = xmin * lsSig;
    const double xplus_s = xplus * lsSig;

    if     (xmin_s>0.0 && xplus_s>-xmin_s) dx = ddwenox(f,1.0);
    else if(xplus_s<0.0 && xmin_s<-xplus_s) dx = ddwenox(f,-1.0);

// y
    if constexpr (Is3D)
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

    const double dnorm2 = dx*dx + dy*dy + dz*dz;
    const double dnorm  = sqrt(dnorm2);

    const double deltax = Is3D ? (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP])
                               : (1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);

    const double denom2 = lsv2 + dnorm2 * deltax * deltax;
    const double sign   = (denom2 > 0.0) ? lsv / sqrt(denom2) : lsSig;

    return sign * (1.0 - dnorm);
}
