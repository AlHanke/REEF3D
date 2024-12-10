/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

void partres::move_RK2(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, turbulence *pturb)
{
    count_particles(p,a,pgc,s);
    pressure_gradient(p,a,pgc,s);
    
    // RK step 1
    stress_tensor(p,pgc,s);
    stress_gradient(p,a,pgc,s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {        
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
        
        if(p->Q11==2)
        advec_mppic(p, a, P, s, pturb, 
                        P.X, P.Y, P.Z, P.U, P.V, P.W,
                        F, G, H, 1.0);
                                         
        // Velocity update
        P.URK1[n] = P.U[n] + p->dtsed*F;
        P.VRK1[n] = P.V[n] + p->dtsed*G;
        P.WRK1[n] = P.W[n] + p->dtsed*H;
        
        // Position update
        P.XRK1[n] = P.X[n] + p->dtsed*P.URK1[n];
        P.YRK1[n] = P.Y[n] + p->dtsed*P.VRK1[n];
        P.ZRK1[n] = P.Z[n] + p->dtsed*P.WRK1[n];
    }
    
    // parallel transfer
    P.xchange(p,pgc,bedch,1);

    // cellSum update
    cellSum_full_update(p,pgc,1);
    
    bedchange_update(p,pgc,s,1);
    // bedchange(p,a,pgc,s,1);

    boundcheck(p,a,pgc,s,1);

    // RK step 2
    stress_tensor(p, pgc, s);
    stress_gradient(p,a,pgc,s);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        if(p->Q11==1)
        advec_plain(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        if(p->Q11==2)
        advec_mppic(p, a, P, s, pturb, 
                        P.XRK1, P.YRK1, P.ZRK1, P.URK1, P.VRK1, P.WRK1,
                        F, G, H, 0.5);
                        
        // Velocity update
        P.U[n] = 0.5*P.U[n] + 0.5*P.URK1[n] + 0.5*p->dtsed*F;
        P.V[n] = 0.5*P.V[n] + 0.5*P.VRK1[n] + 0.5*p->dtsed*G;
        P.W[n] = 0.5*P.W[n] + 0.5*P.WRK1[n] + 0.5*p->dtsed*H;

        if(P.U[n]!=0 || P.V[n]!=0 || P.W[n]!=0)
        {
            P.Flag[n] = MOVING;
        }
        else
        {
            P.Flag[n] = ACTIVE;
        }
        P.Test2[n] = P.U[n];
        
        // Position update
        P.X[n] = 0.5*P.X[n] + 0.5*P.XRK1[n] + 0.5*p->dtsed*P.U[n];
        P.Y[n] = 0.5*P.Y[n] + 0.5*P.YRK1[n] + 0.5*p->dtsed*P.V[n];
        P.Z[n] = 0.5*P.Z[n] + 0.5*P.ZRK1[n] + 0.5*p->dtsed*P.W[n];

        // Resseed relax boundary
        i=p->posc_i(P.XRK1[n]);
        if(P.X[n]>relax_boundary+p->DXN[IP] && P.XRK1[n]<relax_boundary+p->DXN[IP])
        {
            j = p->posc_j(P.YRK1[n]);
            k = p->posc_k(P.ZRK1[n]);
            // seed_particle(p);
        }
    }
    
    // parallel transfer
    P.xchange(p, pgc,bedch,2);
    
    // cellSum update
    cellSum_full_update(p,pgc,2);
    
    bedchange_update(p,pgc,s,2);
    // bedchange(p,a,pgc,s,2);
    
    

    boundcheck(p,a,pgc,s,2);

    ALOOP
    {
        a->test(i,j,k) = p->W1 * 0.5 * PI/8.0 * pow(P.d50,2)  * (pow((a->u(i,j,k)+a->u(i+1,j,k))/2.0,2.0)+pow((a->v(i,j,k)+a->v(i,j+1,k))/2.0,2.0));
        a->test2(i,j,k) = p->Q30*(p->S22-p->W1)*fabs(p->W22)*PI*pow(P.d50, 3.0)/6.0;
        a->test3(i,j,k) = a->test(i,j,k)-a->test2(i,j,k);
        a->test4(i,j,k) = a->test3(i,j,k)>0?1:0;
    }

    // for(n=0;n<P.index;++n)
    // if(P.Flag[n]>=ACTIVE)
    // if(P.X[n]>p->global_xmax || P.Y[n]> p->global_ymax)
    // P.remove(n);
}