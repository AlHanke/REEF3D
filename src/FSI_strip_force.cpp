/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2026 Tobias Martin

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"FSI_strip.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

void fsi_strip::interpolate_vel(lexer* p, fdm* a, ghostcell* pgc, field& uvel, field& vvel, field& wvel)
{
    int ii, jj, kk;
    double dx, dy, dz, dist, D;

    for (int eI = 0; eI < Ne; eI++)
    {
        lagrangeVel[eI] = Eigen::MatrixXd::Zero(3,lagrangePoints[eI].cols());   
    
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            const Eigen::Vector3d& coordI = lagrangePoints[eI].col(pI);

            if 
            (
                coordI(0) >= xstart[p->mpirank] && coordI(0) < xend[p->mpirank] &&
                coordI(1) >= ystart[p->mpirank] && coordI(1) < yend[p->mpirank] &&
                coordI(2) >= zstart[p->mpirank] && coordI(2) < zend[p->mpirank]
            )
            {
                ii = p->posc_i(coordI(0));
                jj = p->posc_j(coordI(1));
                kk = p->posc_k(coordI(2));
                
                dx = p->DXN[IIP];
                dy = p->DYN[JJP];
                dz = p->DZN[KKP];

                for (int i = ii - 2; i <= ii + 2; i++)
                {
                    for (int j = jj - 2; j <= jj + 2; j++)
                    {
                        for (int k = kk - 2; k <= kk + 2; k++)
                        {
                            dist = (p->XN[IP1] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[JP] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[KP] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                            
                            lagrangeVel[eI](0,pI) += uvel(i,j,k)*D;

                            dist = (p->XP[IP] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YN[JP1] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZP[KP] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                                
                            lagrangeVel[eI](1,pI) += vvel(i,j,k)*D;
                            
                            dist = (p->XP[IP] - coordI(0))/dx;
                            D = kernel_roma(dist);
                            dist = (p->YP[JP] - coordI(1))/dy;
                            D *= kernel_roma(dist);
                            dist = (p->ZN[KP1] - coordI(2))/dz;
                            D *= kernel_roma(dist);
                             
                            lagrangeVel[eI](2,pI) += wvel(i,j,k)*D;
                        }
                    }
                }
            }
        }
    }
    
    starttime=pgc->timer();
    for (int eI = 0; eI < Ne; eI++)
    {
        for (int pI = 0; pI < lagrangePoints[eI].cols(); pI++)
        {
            lagrangeVel[eI].col(pI) << pgc->globalsum(lagrangeVel[eI](0,pI)), pgc->globalsum(lagrangeVel[eI](1,pI)), pgc->globalsum(lagrangeVel[eI](2,pI));
        }
    }
    
    if(p->mpirank==0)
    cout<<"FSI_sync time: "<<pgc->timer()-starttime<<endl;
}