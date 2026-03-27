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

#ifndef FLUX_HJ_CDS2_H_
#define FLUX_HJ_CDS2_H_

#include"flux.h"
#include"increment.h"

#include"field.h"
#include"fdm.h"

class flux_HJ_CDS2 : public flux, public increment
{
public:
    flux_HJ_CDS2() = default;
    virtual ~flux_HJ_CDS2() = default;

    inline void u_flux(fdm*, int ipol, const field& uvel, double &uflux1, double&) final
    { u_flux_impl(ipol, uvel, uflux1); }

    inline void v_flux(fdm*, int ipol, const field& vvel, double &vflux1, double&) final
    { v_flux_impl(ipol, vvel, vflux1); }

    inline void w_flux(fdm*, int ipol, const field& wvel, double &wflux1, double&) final
    { w_flux_impl(ipol, wvel, wflux1); }

    #if USE_AMREX
    inline void u_flux(fdm*, int ipol, const LocalArr4Const& uvel, double &uflux1, double&) final
    { u_flux_impl(ipol, uvel, uflux1); }

    inline void v_flux(fdm*, int ipol, const LocalArr4Const& vvel, double &vflux1, double&) final
    { v_flux_impl(ipol, vvel, vflux1); }

    inline void w_flux(fdm*, int ipol, const LocalArr4Const& wvel, double &wflux1, double&) final
    { w_flux_impl(ipol, wvel, wflux1); }
    #endif

private:
    template<typename GenericField>
    inline void u_flux_impl(int ipol, GenericField& uvel, double &uflux1)
    {
        if(ipol==1)
        {
            uflux1 = uvel(i,j,k);
        }
        else if(ipol==2)
        {
            uflux1 = 0.25*(uvel(i,j,k) + uvel(i,j+1,k) + uvel(i-1,j,k) + uvel(i-1,j+1,k));
        }
        else if(ipol==3)
        {
            uflux1 = 0.25*(uvel(i,j,k) + uvel(i,j,k+1) + uvel(i-1,j,k) + uvel(i-1,j,k+1));
        }
        else if(ipol==4)
        {
            uflux1 = 0.5*(uvel(i,j,k) + uvel(i-1,j,k));
        }
    }

    template<typename GenericField>
    inline void v_flux_impl(int ipol, GenericField& vvel, double &vflux1)
    {
        if(ipol==1)
        {
            vflux1 = 0.25*(vvel(i,j,k) + vvel(i+1,j,k) + vvel(i,j-1,k) + vvel(i+1,j-1,k));
        }
        else if(ipol==2)
        {
            vflux1 = vvel(i,j,k);
        }
        else if(ipol==3)
        {
            vflux1 = 0.25*(vvel(i,j,k) + vvel(i,j,k+1) + vvel(i,j-1,k) + vvel(i,j-1,k+1));
        }
        else if(ipol==4)
        {
            vflux1 = 0.5*(vvel(i,j,k) + vvel(i,j-1,k));
        }
    }

    template<typename GenericField>
    inline void w_flux_impl(int ipol, GenericField& wvel, double &wflux1)
    {
        if(ipol==1)
        {
            wflux1 = 0.25*(wvel(i,j,k) + wvel(i+1,j,k) + wvel(i+1,j,k-1) + wvel(i,j,k-1));
        }
        else if(ipol==2)
        {
            wflux1 = 0.25*(wvel(i,j,k) + wvel(i,j+1,k) + wvel(i,j+1,k-1) + wvel(i,j,k-1));
        }
        else if(ipol==3)
        {
            wflux1 = wvel(i,j,k);
        }
        else if(ipol==4)
        {
            wflux1 = 0.5*(wvel(i,j,k) + wvel(i,j,k-1));
        }
    }
};

#endif
