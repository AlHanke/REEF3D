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

#ifndef FLUX_FACE_FOU_VRANS_H_
#define FLUX_FACE_FOU_VRANS_H_

#include"flux.h"
#include"increment.h"

#include"field.h"
#include"lexer.h"
#include"fdm.h"

class flux_face_FOU_vrans final : public flux, public increment
{
public:
	flux_face_FOU_vrans (lexer *p);
	virtual ~flux_face_FOU_vrans() = default;

    inline void u_flux(fdm* a, int ipol, field& uvel, double &uflux1, double &uflux2) override final
    { u_flux_impl(a, ipol, uvel, uflux1, uflux2); }

    inline void v_flux(fdm* a, int ipol, field& vvel, double &vflux1, double &vflux2) override final
    { v_flux_impl(a, ipol, vvel, vflux1, vflux2); }

    inline void w_flux(fdm* a, int ipol, field& wvel, double &wflux1, double &wflux2) override final
    { w_flux_impl(a, ipol, wvel, wflux1, wflux2); }

    #if USE_AMREX
    inline void u_flux(fdm* a, int ipol, const amrex::Array4<const amrex::Real>& uvel, double &uflux1, double &uflux2) override final
    { u_flux_impl(a, ipol, uvel, uflux1, uflux2); }

    inline void v_flux(fdm* a, int ipol, const amrex::Array4<const amrex::Real>& vvel, double &vflux1, double &vflux2) override final
    { v_flux_impl(a, ipol, vvel, vflux1, vflux2); }

    inline void w_flux(fdm* a, int ipol, const amrex::Array4<const amrex::Real>& wvel, double &wflux1, double &wflux2) override final
    { w_flux_impl(a, ipol, wvel, wflux1, wflux2); }
    #endif

private:
    template<typename GenericField>
    void u_flux_impl(fdm* a, int ipol, GenericField& uvel, double &uflux1, double &uflux2)
    {
        if(ipol==1)
        {
            if(p->flag1[Im1JK]>0)
            {
                if(0.5*(uvel(i,j,k)+uvel(i-1,j,k)) >= 0.0)
                uflux1 = uvel(i-1,j,k)*(1.0/a->porosity(i-1,j,k));
                else
                uflux1 = uvel(i,j,k)*(1.0/a->porosity(i-1,j,k));
            }
            else if(p->flag1[Im1JK]<0)
            uflux1 = 0.5*(uvel(i,j,k)+uvel(i-1,j,k))*(1.0/a->porosity(i-1,j,k));

            if(p->flag1[Ip1JK]>0)
            {
                if(0.5*(uvel(i,j,k)+uvel(i+1,j,k)) >= 0.0)
                uflux2 = uvel(i,j,k)*(1.0/a->porosity(i,j,k));
                else
                uflux2 = uvel(i+1,j,k)*(1.0/a->porosity(i,j,k));
            }
            else if(p->flag1[Ip1JK]<0)
            uflux2 = 0.5*(uvel(i,j,k)+uvel(i+1,j,k))*(1.0/a->porosity(i,j,k));
        }
        else if(ipol==2)
        {
            uflux1 = 0.5*(uvel(i-1,j,k)+uvel(i-1,j+1,k))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j+1,k)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
            uflux2 = 0.5*(uvel(i,j,k)+uvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j,k)+a->porosity(i+1,j+1,k))));
        }
        else if(ipol==3)
        {
            uflux1 = 0.5*(uvel(i-1,j,k)+uvel(i-1,j,k+1))*(1.0/(0.25*(a->porosity(i-1,j,k)+a->porosity(i-1,j,k+1)+a->porosity(i,j,k)+a->porosity(i,j,k+1))));
            uflux2 = 0.5*(uvel(i,j,k)+uvel(i,j,k+1))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k)+a->porosity(i+1,j,k+1))));
        }
        else if(ipol==4)
        {
            uflux1 = uvel(i-1,j,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i-1,j,k))));
            uflux2 = uvel(i,j,k)*(1.0/(0.5*(a->porosity(i+1,j,k)+a->porosity(i,j,k))));
        }
    }

    template<typename GenericField>
    void v_flux_impl(fdm* a, int ipol, GenericField& vvel, double &vflux1, double &vflux2)
    {
        if(ipol==1)
        {
            vflux1 = 0.5*(vvel(i,j-1,k)+vvel(i+1,j-1,k))*(1.0/(0.25*(a->porosity(i,j-1,k)+a->porosity(i+1,j-1,k)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
            vflux2 = 0.5*(vvel(i,j,k)+vvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j+1,k)+a->porosity(i+1,j+1,k))));
        }
        else if(ipol==2)
        {
            if(p->flag2[IJm1K]>0)
            {
                if(0.5*(vvel(i,j,k)+vvel(i,j-1,k)) >= 0.0)
                vflux1 = vvel(i,j-1,k)*(1.0/a->porosity(i,j-1,k));
                else
                vflux1 = vvel(i,j,k)*(1.0/a->porosity(i,j-1,k));
            }
            else if(p->flag2[IJm1K]<0)
            vflux1 = 0.5*(vvel(i,j,k)+vvel(i,j-1,k))*(1.0/a->porosity(i,j-1,k));

            if(p->flag2[IJp1K]>0)
            {
                if(0.5*(vvel(i,j,k)+vvel(i,j+1,k)) >= 0.0)
                vflux2 = vvel(i,j,k)*(1.0/a->porosity(i,j,k));
                else
                vflux2 = vvel(i,j+1,k)*(1.0/a->porosity(i,j,k));
            }
            else if(p->flag2[IJp1K]<0)
            vflux2 = 0.5*(vvel(i,j,k)+vvel(i,j+1,k))*(1.0/a->porosity(i,j,k));
        }
        else if(ipol==3)
        {
            vflux1 = 0.5*(vvel(i,j-1,k)+vvel(i,j-1,k+1))*(1.0/(0.25*(a->porosity(i,j-1,k)+a->porosity(i,j-1,k+1)+a->porosity(i,j,k)+a->porosity(i,j,k+1))));
            vflux2 = 0.5*(vvel(i,j,k)+vvel(i,j,k+1))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j,k+1)+a->porosity(i,j+1,k)+a->porosity(i,j+1,k+1))));
        }
        else if(ipol==4)
        {
            vflux1 = vvel(i,j-1,k)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j-1,k))));
            vflux2 = vvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j+1,k)+a->porosity(i,j,k))));
        }
    }

    template<typename GenericField>
    void w_flux_impl(fdm* a, int ipol, GenericField& wvel, double &wflux1, double &wflux2)
    {
        if(ipol==1)
        {
            wflux1 = 0.5*(wvel(i,j,k-1)+wvel(i+1,j,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i+1,j,k-1)+a->porosity(i,j,k)+a->porosity(i+1,j,k))));
            wflux2 = 0.5*(wvel(i,j,k)+wvel(i+1,j,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i+1,j,k)+a->porosity(i,j,k+1)+a->porosity(i+1,j,k+1))));
        }
        else if(ipol==2)
        {
            wflux1 = 0.5*(wvel(i,j,k-1)+wvel(i,j+1,k-1))*(1.0/(0.25*(a->porosity(i,j,k-1)+a->porosity(i,j+1,k-1)+a->porosity(i,j,k)+a->porosity(i,j+1,k))));
            wflux2 = 0.5*(wvel(i,j,k)+wvel(i,j+1,k))*(1.0/(0.25*(a->porosity(i,j,k)+a->porosity(i,j+1,k)+a->porosity(i,j,k+1)+a->porosity(i,j+1,k+1))));
        }
        else if(ipol==3)
        {
            if(p->flag3[IJKm1]>0)
            {
                if(0.5*(wvel(i,j,k)+wvel(i,j,k-1))>=0.0)
                wflux1 = wvel(i,j,k-1)*(1.0/a->porosity(i,j,k-1));
                else
                wflux1 = wvel(i,j,k)*(1.0/a->porosity(i,j,k-1));
            }
            else if(p->flag3[IJKm1]<0)
            wflux1 = 0.5*(wvel(i,j,k)+wvel(i,j,k-1))*(1.0/a->porosity(i,j,k-1));

            if(p->flag3[IJKp1]>0)
            {
                if(0.5*(wvel(i,j,k)+wvel(i,j,k+1))>=0.0)
                wflux2 = wvel(i,j,k)*(1.0/a->porosity(i,j,k));
                else
                wflux2 = wvel(i,j,k+1)*(1.0/a->porosity(i,j,k));
            }
            else if(p->flag3[IJKp1]<0)
            wflux2 = 0.5*(wvel(i,j,k)+wvel(i,j,k+1))*(1.0/a->porosity(i,j,k));
        }
        else if(ipol==4)
        {
            wflux1 = wvel(i,j,k-1)*(1.0/(0.5*(a->porosity(i,j,k)+a->porosity(i,j,k-1))));
            wflux2 = wvel(i,j,k)*(1.0/(0.5*(a->porosity(i,j,k+1)+a->porosity(i,j,k))));
        }
    }

    lexer *p;
};

#endif
