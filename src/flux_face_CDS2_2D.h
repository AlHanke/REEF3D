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

#ifndef FLUX_FACE_CDS2_2D_H_
#define FLUX_FACE_CDS2_2D_H_

#include"flux.h"
#include"increment.h"

#include"field.h"

class flux_face_CDS2_2D final : public flux, public increment
{
public:
    flux_face_CDS2_2D() = default;
    virtual ~flux_face_CDS2_2D() = default;

    inline void u_flux(fdm*, int ipol, const field& uvel, double &uflux1, double &uflux2) override final
    { u_flux_impl(ipol, uvel, uflux1, uflux2); }

    inline void v_flux(fdm*, int, const field&, double &vflux1, double &vflux2) override final
    { vflux1 = 0.0; vflux2 = 0.0; }

    inline void w_flux(fdm*, int ipol, const field& wvel, double &wflux1, double &wflux2) override final
    { w_flux_impl(ipol, wvel, wflux1, wflux2); }

    #if USE_AMREX
    inline void u_flux(fdm*, int ipol, const LocalArr4Const& uvel, double &uflux1, double &uflux2) override final
    { u_flux_impl(ipol, uvel, uflux1, uflux2); }

    inline void v_flux(fdm*,int, const LocalArr4Const&, double &vflux1, double &vflux2) override final
    { vflux1 = 0.0; vflux2 = 0.0; }

    inline void w_flux(fdm*, int ipol, const LocalArr4Const& wvel, double &wflux1, double &wflux2) override final
    { w_flux_impl(ipol, wvel, wflux1, wflux2); }
    #endif

private:
    template<typename GenericField>
    inline void u_flux_impl(int ipol, GenericField& uvel, double &uflux1, double &uflux2)
    {
        if(ipol==1)
        {
            uflux1 = 0.5*(uvel(i,j,k)+uvel(i-1,j,k));
            uflux2 = 0.5*(uvel(i,j,k)+uvel(i+1,j,k));
        }
        else if(ipol==2)
        {
            uflux1 = uvel(i-1,j,k);
            uflux2 = uvel(i,j,k);
        }
        else if(ipol==3)
        {
            uflux1 = 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1));
            uflux2 = 0.5*(uvel(i,j,k)+uvel(i,j,k+1));
        }
        else if(ipol==4)
        {
            uflux1 = uvel(i-1,j,k);
            uflux2 = uvel(i,j,k);
        }
    }
    template<typename GenericField>
    inline void w_flux_impl(int ipol, GenericField& wvel, double &wflux1, double &wflux2)
    {
        if(ipol==1)
        {
            wflux1 = 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1));
            wflux2 = 0.5*(wvel(i,j,k)+wvel(i+1,j,k));
        }
        else if(ipol==2)
        {
            wflux1 = wvel(i,j,k-1);
            wflux2 = wvel(i,j,k);
        }
        else if(ipol==3)
        {
            wflux1 = 0.5*(wvel(i,j,k)+wvel(i,j,k-1));
            wflux2 = 0.5*(wvel(i,j,k)+wvel(i,j,k+1));
        }
        else if(ipol==4)
        {
            wflux1 = wvel(i,j,k-1);
            wflux2 = wvel(i,j,k);
        }
    }
};

#endif
