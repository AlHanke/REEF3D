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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

#include <map>

void partres::cellSum_update(lexer *p, ghostcell *pgc, sediment_fdm *s, int mode)
{
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        // step 1
        if(mode==1)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        cellSum(i,j,k) -= P.ParcelFactor;
        bedch(i,j) -= P.ParcelFactor;
        
        
        // step 2
        if(mode==1)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        cellSum(i,j,k) += P.ParcelFactor;
        bedch(i,j) += P.ParcelFactor;
    }
    
    pgc->gcsl_start4(p,bedch,1);
}

void partres::cellSum_full_update(lexer *p, ghostcell *pgc, int mode)
{
    ALOOP
    cellSum(i,j,k) = 0.0;
    pgc->start4a(p,cellSum,1);
    
    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]>=ACTIVE)
    {
        if(mode==1)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        
        Sx = (p->XN[IP1] - P.XRK1[n])/(p->XN[IP1] - p->XN[IP]);
        Sy = (p->YN[JP1] - P.YRK1[n])/(p->YN[JP1] - p->YN[JP]);
        Sz = (p->ZN[KP1] - P.ZRK1[n])/(p->ZN[KP1] - p->ZN[KP]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        
        Sx = (p->XN[IP1] - P.X[n])/(p->XN[IP1] - p->XN[IP]);
        Sy = (p->YN[JP1] - P.Y[n])/(p->YN[JP1] - p->YN[JP]);
        Sz = (p->ZN[KP1] - P.Z[n])/(p->ZN[KP1] - p->ZN[KP]);
        }
        
        cellSum(i,j,k) += P.ParcelFactor * Sx*Sy*Sz;
        cellSum(i+1,j,k) += P.ParcelFactor * (1.0-Sx)*Sy*Sz;
        cellSum(i+1,j+1,k) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*Sz;
        cellSum(i,j+1,k) += P.ParcelFactor * Sx*(1.0-Sy)*Sz;
        cellSum(i,j,k+1) += P.ParcelFactor * Sx*Sy*(1.0-Sz);
        cellSum(i+1,j,k+1) += P.ParcelFactor * (1.0-Sx)*Sy*(1.0-Sz);
        cellSum(i+1,j+1,k+1) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*(1.0-Sz);
        cellSum(i,j+1,k+1) += P.ParcelFactor * Sx*(1.0-Sy)*(1.0-Sz);
    }
    pgc->start4a_sum(p,cellSum,1);
    if(p->nb5==-2)
    {
        k=0;
        ILOOP
        JLOOP
        cellSum(i,j,k) *= 2;
    }
    if(p->nb1==-2)
    {
        i=0;
        KLOOP
        JLOOP
        cellSum(i,j,k) *= 2;
    }
    if(p->nb3==-2)
    {
        j=0;
        ILOOP
        KLOOP
        cellSum(i,j,k) *= 2;
    }
    if(p->nb1!=-2 && p->nb3!=-2)
    {
        if(p->nb5==-2)
            cellSum(0,0,0) *= 4.0/3.0;
    }
    pgc->start4a(p,cellSum,1);
}

void partres::cellSum_gradient_redistirbution(lexer *p, ghostcell *pgc, int mode)
{
    for(int iter=0; iter<5; iter++)
    {
        cellSum_full_update(p,pgc,mode);
        std::map<std::tuple<int,int,int>,std::vector<int>> cellSum_map;

        for(size_t n=0;n<P.index;n++)
        if(P.Flag[n]>=ACTIVE)
        {
            if(mode==1)
            {
                i=p->posc_i(P.XRK1[n]);
                j=p->posc_j(P.YRK1[n]);
                k=p->posc_k(P.ZRK1[n]);
            }
            
            if(mode==2)
            {
                i=p->posc_i(P.X[n]);
                j=p->posc_j(P.Y[n]);
                k=p->posc_k(P.Z[n]);
            }

            cellSum_map[std::make_tuple(i,j,k)].push_back(n);
        }
        ALOOP
            if(Ts(i,j,k)>Tc)
            {
                double grad = cellSum(i,j,k)-cellSum(i+1,j,k);
                int dir = 4;
                if(cellSum(i,j,k)-cellSum(i-1,j,k)>grad)
                {
                    grad = cellSum(i,j,k)-cellSum(i-1,j,k);
                    dir = 1;
                }
                if(cellSum(i,j,k)-cellSum(i,j+1,k)>grad)
                {
                    grad = cellSum(i,j,k)-cellSum(i,j+1,k);
                    dir = 2;
                }
                if(cellSum(i,j,k)-cellSum(i,j-1,k)>grad)
                {
                    grad = cellSum(i,j,k)-cellSum(i,j-1,k);
                    dir = 3;
                }
                if(cellSum(i,j,k)-cellSum(i,j,k+1)>grad)
                {
                    grad = cellSum(i,j,k)-cellSum(i,j,k+1);
                    dir = 6;
                }
                if(cellSum(i,j,k)-cellSum(i,j,k-1)>grad)
                {
                    grad = cellSum(i,j,k)-cellSum(i,j,k-1);
                    dir = 5;
                }
                if(0<grad)
                {
                    size_t moveIndex;
                    size_t cellIndex=-1;
                    double maxDist;
                    double pos;
                    if (dir==4)
                    {
                        maxDist = p->DXP[IP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.XRK1[n];
                            if(mode==2)
                                pos = P.X[n];
                            if(p->XN[IP1]-pos<maxDist)
                            {
                                maxDist = p->XN[IP1]-pos;
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.XRK1[moveIndex] += p->DXP[IP];
                            P.X[moveIndex] += p->DXP[IP];
                            cellSum(i+1,j,k) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i+1,j,k)].push_back(moveIndex);
                        }
                    }
                    if (dir==1)
                    {
                        maxDist = p->DXP[IP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.XRK1[n];
                            if(mode==2)
                                pos = P.X[n];
                            if(pos-p->XN[IP]<maxDist)
                            {
                                maxDist = pos-p->XN[IP];
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.XRK1[moveIndex] -= p->DXP[IM1];
                            P.X[moveIndex] -= p->DXP[IM1];
                            cellSum(i-1,j,k) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i-1,j,k)].push_back(moveIndex);
                        }
                    }
                    if (dir==2)
                    {
                        maxDist = p->DYP[JP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.YRK1[n];
                            if(mode==2)
                                pos = P.Y[n];
                            if(p->YN[JP1]-pos<maxDist)
                            {
                                maxDist = p->YN[JP1]-pos;
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.YRK1[moveIndex] += p->DYP[JP];
                            P.Y[moveIndex] += p->DYP[JP];
                            cellSum(i,j+1,k) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i,j+1,k)].push_back(moveIndex);
                        }
                    }
                    if (dir==3)
                    {
                        maxDist = p->DYP[JP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.YRK1[n];
                            if(mode==2)
                                pos = P.Y[n];
                            if(pos-p->YN[JP]<maxDist)
                            {
                                maxDist = pos-p->YN[JP];
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.YRK1[moveIndex] -= p->DYP[JM1];
                            P.Y[moveIndex] -= p->DYP[JM1];
                            cellSum(i,j-1,k) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i,j-1,k)].push_back(moveIndex);
                        }
                    }
                    if (dir==6)
                    {
                        maxDist = p->DZP[KP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.ZRK1[n];
                            if(mode==2)
                                pos = P.Z[n];
                            if(p->ZN[KP1]-pos<maxDist)
                            {
                                maxDist = pos-p->ZN[KP1];
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.ZRK1[moveIndex] += p->DZP[KP];
                            P.Z[moveIndex] += p->DZP[KP];
                            cellSum(i,j,k+1) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i,j,k+1)].push_back(moveIndex);
                        }
                    }
                    if (dir==5)
                    {
                        maxDist = p->DZP[KP];
                        for(size_t n=0;n<cellSum_map[std::make_tuple(i,j,k)].size();n++)
                        {
                            if(mode==1)
                                pos = P.ZRK1[n];
                            if(mode==2)
                                pos = P.Z[n];
                            if(pos-p->ZN[KP]<maxDist)
                            {
                                maxDist = pos-p->ZN[KP];
                                moveIndex = cellSum_map[std::make_tuple(i,j,k)][n];
                                cellIndex = n;
                            }
                        }
                        if(cellIndex!=-1)
                        {
                            P.ZRK1[moveIndex] -= p->DZP[KM1];
                            P.Z[moveIndex] -= p->DZP[KM1];
                            cellSum(i,j,k-1) += P.ParcelFactor;
                            cellSum(i,j,k) -= P.ParcelFactor;
                            cellSum_map[std::make_tuple(i,j,k)].erase(cellSum_map[std::make_tuple(i,j,k)].begin()+cellIndex);
                            cellSum_map[std::make_tuple(i,j,k-1)].push_back(moveIndex);
                        }
                    }
                }
            }
        P.xchange(p,pgc,bedch,mode);
    }
    cellSum_full_update(p,pgc,mode);
}
