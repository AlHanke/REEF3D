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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

void initialize::hydrostatic(lexer *p, fdm *a, ghostcell *pgc)
{
    double maxh=0.0;
    BASELOOP
    maxh=MAX(maxh, p->pos_z());

    maxh=pgc->globalmax(maxh);
	
    if(p->F30>0)
    maxh=p->phimean;
	
	if(p->I12==1 && (p->I30==0||p->B90==0))
    BASELOOP
    a->press(i,j,k) = (p->phimean-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

	if(p->I12==2 && (p->I30==0||p->B90==0))
    BASELOOP
    a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22);
	
	if(p->I12==3 && (p->I30==0||p->B90==0))
    BASELOOP
    a->press(i,j,k) = (maxh-p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

    if(p->I12==4 && p->Y9==1)
    {
        #if USE_AMREX
        const double psi = p->psi;
        auto Hface = [&](double phiv)->double {
            if(phiv> psi) return 1.0;
            if(phiv<-psi) return 0.0;
            return 0.5*(1.0 + phiv/psi + (1.0/PI)*sin((PI*phiv)/psi));
        };

        // LEVEL_LOOP
        // // TILE_LOOP
        // ILOOP
        // JLOOP
        // {
        //     // anchor the bottom cell on the analytic value (absolute level is irrelevant to velcorr)
        //     k=0;
        //     a->press0(i,j,k) = a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22);

        //     // integrate the DISCRETE balance upward: press(k+1) = press(k) + W22*DZP*roface
        //     for(k=0; k<KMAX_LOOP; ++k)
        //     {
        //         const double phiface = 0.5*(a->phi(i,j,k) + a->phi(i,j,k+1));
        //         const double H       = Hface(phiface);
        //         const double roface  = p->W1*H + p->W3*(1.0-H);   // matches density_f::roface
        //         a->press0(i,j,k+1) = a->press(i,j,k+1)    = a->press(i,j,k) + p->W22*p->DZP[KP]*roface;
        //     }
        // }

        BASELOOP
        {
            const double dz  = p->amrex_geometry[p->level].CellSize(2);
            const double zlo = p->amrex_geometry[p->level].ProbLo(2);
            const int    gk  = k + p->amr_tile_lo.z;     // GLOBAL k on this level

            const double zc   = p->pos_z();              // this cell's centre z
            const double phic = a->phi(i,j,k);           // signed distance, dphi/dz = -1

            // reference = column bottom cell (gk=0): deep water, roface==W1, analytic==discrete
            const double zc0  = zlo + 0.5*dz;
            const double phi0 = phic + (zc - zc0);       // phi at bottom via unit gradient
            double press = phi0*p->W1*fabs(p->W22) + p->I55;

            // discrete balance, face by face, up to this cell -- all evaluated ANALYTICALLY
            // (no field reads off-tile), so it is order-independent across tiles and levels
            for(int m=0; m<gk; ++m)
            {
                const double zf  = zlo + double(m+1)*dz;     // face between cells m and m+1
                const double phif= phic + (zc - zf);         // phi at that face, unit gradient
                const double H   = Hface(phif);
                const double rof = p->W1*H + p->W3*(1.0-H);  // == density_f::roface there
                press += p->W22*dz*rof;                      // W22 < 0
            }

            a->press0(i,j,k) = a->press(i,j,k) = press;
        }

        pgc->start4(p,a->press,40);   // fill halos / C-F ghosts on every level
        pgc->start4(p,a->press0,40);
        #endif
    }

    if(p->I56==1)
    BASELOOP
    {
    if(a->phi(i,j,k)<0.0)
    a->press(i,j,k)=0.0;
    }
	
    BASELOOP
    a->press(i,j,k)+=p->I55;
	
	if(p->I12==2 && p->I30==0)
    GC4LOOP
	{
	i = p->gcb4[n][0];
	j = p->gcb4[n][1];
	k = p->gcb4[n][2];
	
    a->press(i,j,k) = a->phi(i,j,k)*a->ro(i,j,k)*fabs(p->W22) + p->I55;
	}

    
    
    
	

}





