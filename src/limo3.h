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

#ifndef LIMO3_H_
#define LIMO3_H_

#include"fluxlim.h"
#include"increment.h"
#include"lexer.h"

class limo3 final : public fluxlim, public increment
{
public:
    limo3(lexer *p) : delta(p->DXM), radius (0.1), eps(1.0e-9*p->DXM) {};
    virtual ~limo3() = default;

protected:
    inline double phi_impl(double vn1, double vn2, double vq1, double vq2, double) override final
    {
        double d1 = vn1 - vn2;
        double d2 = vq1 - vq2;

        double r = d1 / (fabs(d2)>1.0e-10 ? d2 : 1.0e20);
        double eta = (d1*d1 + d2*d2) / pow(radius*delta, 2.0);
        double phihat = max(0.0, min((2.0+r)/3.0, max(-0.5*r, min(2.0*r, (2.0+r)/3.0, 1.6))));

        if(eta <= 1.0-eps)           return (2.0+r)/3.0;
        else if(eta >= 1.0+eps)      return phihat;
        else return 0.5*((1.0-(eta-1.0)/eps)*(2.0+r)/3.0 + (1.0+(eta-1.0)/eps)*phihat);
    };

private:
    inline double max(double val1, double val2, double val3)
    {
        double maxi;

        maxi=val1;

        if(maxi<val2)
        maxi=val2;

        if(maxi<val3)
        maxi=val3;

        if(maxi<0.0)
        maxi=0.0;

        return maxi;
    };
    inline double max(double val1, double val2)
    {
        double maxi;

        maxi=val1;

        if(maxi<val2)
        maxi=val2;

        if(maxi<0.0)
        maxi=0.0;

        return maxi;
    };
    inline double min(double val1, double val2, double val3)
    {
        double mini;

        mini=val1;

        if(mini>val2)
        mini=val2;

        if(mini>val3)
        mini=val3;

        if(mini<0.0)
        mini=0.0;

        return mini;
    };
    inline double min(double val1,double val2)
    {
        double mini;

        mini=val1;

        if(mini>val2)
        mini=val2;

        if(mini<0.0)
        mini=0.0;

        return mini;
    };

    const double delta,radius,eps;
};

#endif

