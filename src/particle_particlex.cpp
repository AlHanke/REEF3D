/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<array>

void particle_f::particlex(lexer* p, fdm* a, ghostcell* pgc)
{
    xchange=0;

    size_t maxcap=ceil(p->Q25*double(gpartnum));
    tracers_obj seedling1(maxcap),seedling2(maxcap),seedling3(maxcap),seedling4(maxcap),seedling5(maxcap),seedling6(maxcap);
    tracers_obj Send[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};
    tracers_obj Recv[6]={seedling1,seedling2,seedling3,seedling4,seedling5,seedling6};

    PARTLOOP
        if(PP.Flag[n]==1)
        {
            i = p->posc_i(PP.X[n]);
            j = p->posc_j(PP.Y[n]);
            k = p->posc_k(PP.Z[n]);


            if(p->flag5[IJK]<0 && p->flag5[IJK]>-10)
            {
                switch (p->flag5[IJK])
                {
                    case -1:
                    {
                        Send[0].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -2:
                    {
                        Send[1].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -3:
                    {
                        Send[2].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -4:
                    {
                        Send[3].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -5:
                    {
                        Send[4].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -6:
                    {
                        Send[5].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }
                }
                PP.erase(n);
                ++xchange;
            }
        }

    pgc->para_tracersobj(p,Send,Recv);


    size_t sum=PP.size;
    for(int n=0;n<6;n++)
        sum += Recv[n].size;
    if(sum>PP.capacity)
        PP.reserve(sum);

    for(int n=0;n<6;n++)
    {
        PP.add_obj(&Recv[n]);
        // Send[n].erase_all();
        // Recv[n].erase_all();
    }
}


