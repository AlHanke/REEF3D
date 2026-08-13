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

#include"iowave.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_active_wavegen(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double *UH, double *VH, double *WH, slice &WL)
{
    double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
        
        double etaval=0.0;
        
        // wave theory
        if(p->B92<20 || p->B92>29)
        GCSLIN
        {
        i=p->gcslin[p->level][n].i;
        j=p->gcslin[p->level][n].j;
        
            if(p->A515==1)
            etaval = d->eta(i,j);
            
            if(p->A515==2)
            {
            etaval = d->eta(i,j);
            /*
            WL(i-1,j) = WL(i,j);
            WL(i-2,j) = WL(i,j);
            WL(i-3,j) = WL(i,j);*/
            }
            
            if(p->A515==3)
            {
            etaval = eta(i-1,j);
            
            WL(i-1,j) = etaval + d->depth(i,j);
            WL(i-2,j) = etaval + d->depth(i,j);
            WL(i-3,j) = etaval + d->depth(i,j);
            }
        
        d->eta(i-1,j) = etaval;
        d->eta(i-2,j) = etaval;
        d->eta(i-3,j) = etaval;
        }
        
        // wave maker
        if(p->B92>=20 && p->B92<=29)
        GCSLIN
        {
        i=p->gcslin[p->level][n].i;
        j=p->gcslin[p->level][n].j;
        
            if(p->A515==1)
            etaval = d->eta(i,j);
            
            if(p->A515==2)
            {
            etaval = d->eta(i,j);
            /*
            WL(i-1,j) = WL(i,j);
            WL(i-2,j) = WL(i,j);
            WL(i-3,j) = WL(i,j);*/
            }
            
            if(p->A515==3)
            {
            etaval = eta(i-1,j);
            
            WL(i-1,j) = etaval + d->depth(i,j);
            WL(i-2,j) = etaval + d->depth(i,j);
            WL(i-3,j) = etaval + d->depth(i,j);
            }
        d->eta(i-1,j) = etaval;
        d->eta(i-2,j) = etaval;
        d->eta(i-3,j) = etaval;
        }
        
        
        // wavegen
        count=0;
        // LEVEL_LOOP
        GCINLOOP
        {
        i=p->gcin[p->level][n].i;
        j=p->gcin[p->level][n].j;
        k=p->gcin[p->level][n].k;

        GCIN_TILE(n);

            WETDRYDEEP
            {
            uvel=uval[count]*ramp(p);
            vvel=vval[count]*ramp(p);
            wvel=wval[count]*ramp(p);

                U[Im1JK]=uvel+p->Ui;
                U[Im2JK]=uvel+p->Ui;
                U[Im3JK]=uvel+p->Ui;
                    
                V[Im1JK]=vvel;
                V[Im2JK]=vvel;
                V[Im3JK]=vvel;
                    
                W[Im1JK]=wvel;
                W[Im2JK]=wvel;
                W[Im3JK]=wvel;
            
            uvel=UHval[count]*ramp(p);
            vvel=VHval[count]*ramp(p);
            wvel=WHval[count]*ramp(p);
            
                UH[Im1JK]=uvel;
                UH[Im2JK]=uvel;
                UH[Im3JK]=uvel;
                
                VH[Im1JK]=vvel;
                VH[Im2JK]=vvel;
                VH[Im3JK]=vvel;
                
                WH[Im1JK]=wvel;
                WH[Im2JK]=wvel;
                WH[Im3JK]=wvel;
            
                
                // fsf deviation
                eta_T = wave_eta(p,pgc,x,0.0);
                eta_M = d->eta(i,j); 
                eta_R = eta_T-eta_M;
				
                Uc=eta_R*sqrt(9.81/p->wd);
                
                U[Im1JK] += Uc;
                U[Im2JK] += Uc;
                U[Im3JK] += Uc;
                
                UH[Im1JK] += (eta_R+d->depth(i,j))*Uc;
                UH[Im2JK] += (eta_R+d->depth(i,j))*Uc;
                UH[Im3JK] += (eta_R+d->depth(i,j))*Uc;
         }
         ++count;
		}
        GC_TILE_RESET;

        /*
         if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		{
            // LEVEL_LOOP
            GCINLOOP
            {
            GCIN_TILE(n);

            for(int q=0;q<4;++q)
            {
            i=p->gcin[p->level][n].i+q;
            j=p->gcin[p->level][n].j;
            k=p->gcin[p->level][n].k;

            d->EV[IJK]=MIN(d->EV[IJK],1.0e-4);
            }
            }
            GC_TILE_RESET;
         pgc->start24V(p,d->EV,24);
		}*/
        
}