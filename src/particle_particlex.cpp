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
    int q=0;
    xchange=0;

    for(n=0;n<6;++n)
    {
        pxs[n]=0;
        pxr[n]=0;
    }
    size_t maxcap=ceil(p->Q25*double(gpartnum));
    tracers_obj seedling(maxcap);
    tracers_obj Send[6]={seedling,seedling,seedling,seedling,seedling,seedling};
    tracers_obj Recv[6]={seedling,seedling,seedling,seedling,seedling,seedling};

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
                        posxs[0][pxs[0]+0]=PP.X[n];
                        posxs[0][pxs[0]+1]=PP.Y[n];
                        posxs[0][pxs[0]+2]=PP.Z[n];
                        pxs[0]+=PP.entries;
                        Send[0].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -2:
                    {
                        posxs[1][pxs[1]+0]=PP.X[n];
                        posxs[1][pxs[1]+1]=PP.Y[n];
                        posxs[1][pxs[1]+2]=PP.Z[n];
                        pxs[1]+=PP.entries;
                        Send[1].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -3:
                    {
                        posxs[2][pxs[2]+0]=PP.X[n];
                        posxs[2][pxs[2]+1]=PP.Y[n];
                        posxs[2][pxs[2]+2]=PP.Z[n];
                        pxs[2]+=PP.entries;
                        Send[2].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -4:
                    {
                        posxs[3][pxs[3]+0]=PP.X[n];
                        posxs[3][pxs[3]+1]=PP.Y[n];
                        posxs[3][pxs[3]+2]=PP.Z[n];
                        pxs[3]+=PP.entries;
                        Send[3].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -5:
                    {
                        posxs[4][pxs[4]+0]=PP.X[n];
                        posxs[4][pxs[4]+1]=PP.Y[n];
                        posxs[4][pxs[4]+2]=PP.Z[n];
                        pxs[4]+=PP.entries;
                        Send[4].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }

                    case -6:
                    {
                        posxs[5][pxs[5]+0]=PP.X[n];
                        posxs[5][pxs[5]+1]=PP.Y[n];
                        posxs[5][pxs[5]+2]=PP.Z[n];
                        pxs[5]+=PP.entries;
                        Send[5].add(PP.X[n],PP.Y[n],PP.Z[n],1);
                        break;
                    }
                }
                PP.erase(n);
                ++xchange;
            }
        }
    // cout<<p->mpirank<<"|"<<xchange<<":"<<PP.size<<endl;
    pgc->paratracersobj(p,posxs,posxr,pxs,pxr);
    pgc->paratracersobj2(p,Send,Recv);
    cout<<"done transfer:"<<p->mpirank<<endl;

    for(int n=0;n<6;n++)
    {
        // for(int m=0;m<Recv[n].loopindex;m++)
        //     cout<<p->mpirank<<"|"<<n<<":"<<m<<"|"<<Recv[n].X[m]<<endl;
        PP.add_obj(&Recv[n]);
        Send[n].erase_all();
        Recv[n].erase_all();
    }
    cout<<"done:"<<p->mpirank<<":"<<PP.size<<endl;

    // q=0;
    // while(q<6)
    // {
    //     if(pxr[q]>0)
    //     {
    //         PP.add(posxr[q][pxr[q]-(PP.entries-0)],posxr[q][pxr[q]-(PP.entries-1)],posxr[q][pxr[q]-(PP.entries-2)],1);
    //         pxr[q]-=PP.entries;
    //     }

    //     if(pxr[q]<=0)
    //         ++q;
    // }
    
}


