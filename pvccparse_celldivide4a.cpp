/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"pvccparse.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pvccparse::cell_divide4a(lexer* p, fdm* a, ghostcell *pgc)
{
		// cell 1
        if(pt[0]==0 && pt[4]==0 && cl[3]<=0 && cl[0]<=0 && cl[8]<=0)
        {
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 2
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 3
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 4
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 5
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 6
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 7
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 2
        if(pt[1]==0 && pt[5]==0 && cl[0]<=0 && cl[1]<=0 && cl[9]<=0)
        {
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 3
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 4
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 5
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 6
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 7
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 3
        if(pt[2]==0 && pt[6]==0 && cl[2]<=0 && cl[1]<=0 && cl[10]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[3];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 4
        if(pt[3]==0 && pt[7]==0 && cl[2]<=0 && cl[3]<=0 && cl[11]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 5
        if(pt[0]==0 && pt[4]==0 && cl[11]<=0 && cl[8]<=0 && cl[0]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[5];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 6
        if(pt[1]==0 && pt[5]==0 && cl[8]<=0 && cl[9]<=0 && cl[1]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[0];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 7
        if(pt[2]==0 && pt[6]==0 && cl[10]<=0 && cl[9]<=0 && cl[2]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 8
        if(pt[3]==0 && pt[7]==0 && cl[10]<=0 && cl[11]<=0 && cl[3]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 9
        if(pt[0]==0 && pt[4]==0 && cl[3]<=0 && cl[0]<=0 && cl[11]<=0)
        {
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 2
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 3
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 4
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 5
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 6
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
		
		// tetra 7
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 10
        if(pt[1]==0 && pt[5]==0 && cl[0]<=0 && cl[1]<=0 && cl[8]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[0];
        a->pvccnode[count][3]=nd[3];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[1];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 11
        if(pt[2]==0 && pt[6]==0 && cl[1]<=0 && cl[2]<=0 && cl[9]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 12
        if(pt[3]==0 && pt[7]==0 && cl[2]<=0 && cl[3]<=0 && cl[10]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 13
        if(pt[0]==0 && pt[4]==0 && cl[11]<=0 && cl[8]<=0 && cl[3]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 14
        if(pt[1]==0 && pt[5]==0 && cl[8]<=0 && cl[9]<=0 && cl[0]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[4];
	
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 15
        if(pt[2]==0 && pt[6]==0 && cl[9]<=0 && cl[10]<=0 && cl[1]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[8];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[5];
	
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 16
        if(pt[3]==0 && pt[7]==0 && cl[10]<=0 && cl[11]<=0 && cl[2]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 17
        if(pt[0]==0 && pt[3]==0 && cl[2]<=0 && cl[7]<=0 && cl[4]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[0];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 18
        if(pt[1]==0 && pt[2]==0 && cl[2]<=0 && cl[6]<=0 && cl[0]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[7];
	
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 19
        if(pt[5]==0 && pt[6]==0 && cl[10]<=0 && cl[6]<=0 && cl[5]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 20
        if(pt[4]==0 && pt[7]==0 && cl[7]<=0 && cl[10]<=0 && cl[8]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[6];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 21 
        if(pt[0]==0 && pt[3]==0 && cl[0]<=0 && cl[4]<=0 && cl[7]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[1];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[0];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[1];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 22
        if(pt[1]==0 && pt[2]==0 && cl[0]<=0 && cl[5]<=0 && cl[2]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[6];
	
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 23
        if(pt[5]==0 && pt[6]==0 && cl[8]<=0 && cl[5]<=0 && cl[6]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
	
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 24
        if(pt[4]==0 && pt[7]==0 && cl[4]<=0 && cl[8]<=0 && cl[10]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[6];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 25
        if(pt[0]==0 && pt[3]==0 && cl[2]<=0 && cl[7]<=0 && cl[0]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 26
        if(pt[1]==0 && pt[2]==0 && cl[2]<=0 && cl[6]<=0 && cl[5]<=0)
        {
		
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[6];
	
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[1];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 27
        if(pt[5]==0 && pt[6]==0 && cl[6]<=0 && cl[10]<=0 && cl[8]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
	
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 28
        if(pt[4]==0 && pt[7]==0 && cl[7]<=0 && cl[10]<=0 && cl[4]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		
		// cell 29
        if(pt[0]==0 && pt[3]==0 && cl[4]<=0 && cl[0]<=0 && cl[2]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 30
        if(pt[1]==0 && pt[2]==0 && cl[0]<=0 && cl[5]<=0 && cl[6]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[6];
	
        a->pvccnode[count][4]=nd[1];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
		
        a->pvccnode[count][4]=nd[1];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[4];

        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 31
        if(pt[5]==0 && pt[6]==0 && cl[5]<=0 && cl[8]<=0 && cl[10]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[8];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
	
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 32
        if(pt[4]==0 && pt[7]==0 && cl[4]<=0 && cl[8]<=0 && cl[7]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[1];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[5];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[7];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[5];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 33
        if(pt[0]==0 && pt[1]==0 && cl[4]<=0 && cl[1]<=0 && cl[5]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[3];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[0];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 34
        if(pt[4]==0 && pt[5]==0 && cl[5]<=0 && cl[9]<=0 && cl[11]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[7];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[6];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 35
        if(pt[7]==0 && pt[6]==0 && cl[6]<=0 && cl[9]<=0 && cl[7]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];
	
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[3];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 36
        if(pt[3]==0 && pt[2]==0 && cl[1]<=0 && cl[6]<=0 && cl[3]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[1];
	
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 37
        if(pt[0]==0 && pt[1]==0 && cl[3]<=0 && cl[4]<=0 && cl[5]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 38
        if(pt[4]==0 && pt[5]==0 && cl[4]<=0 && cl[9]<=0 && cl[11]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=nd[3];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[8];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 39
        if(pt[7]==0 && pt[6]==0 && cl[7]<=0 && cl[11]<=0 && cl[6]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];
	
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[3];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 40
        if(pt[3]==0 && pt[2]==0 && cl[3]<=0 && cl[7]<=0 && cl[1]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[1];
	
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[0];
		
        a->pvccnode[count][4]=nd[3];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[3];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 41
        if(pt[0]==0 && pt[1]==0 && cl[1]<=0 && cl[5]<=0 && cl[3]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[5];
        a->pvccnode[count][3]=nd[3];
	
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[5];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 42
        if(pt[4]==0 && pt[5]==0 && cl[5]<=0 && cl[9]<=0 && cl[4]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[6];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 43
        if(pt[7]==0 && pt[6]==0 && cl[6]<=0 && cl[9]<=0 && cl[11]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[2];
        a->pvccnode[count][5]=nd[5];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];
	
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[6];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[4];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 44
        if(pt[3]==0 && pt[2]==0 && cl[1]<=0 && cl[6]<=0 && cl[7]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[0];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[7];
	
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[3];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[2];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[2];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 45
        if(pt[0]==0 && pt[1]==0 && cl[3]<=0 && cl[4]<=0 && cl[1]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[6];
        a->pvccnode[count][4]=nd[7];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[6];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[2];
	
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[2];
        a->pvccnode[count][2]=nd[3];
        a->pvccnode[count][3]=nd[4];
		
        a->pvccnode[count][4]=nd[0];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[0];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 46
        if(pt[4]==0 && pt[5]==0 && cl[4]<=0 && cl[11]<=0 && cl[5]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[0];
	
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[5];
        a->pvccnode[count][3]=nd[1];
		
        a->pvccnode[count][4]=nd[8];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[6];


        a->pvccnode[count][3]=nd[8];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 47
        if(pt[7]==0 && pt[6]==0 && cl[7]<=0 && cl[11]<=0 && cl[9]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[4];
        a->pvccnode[count][2]=nd[3];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[5];
        a->pvccnode[count][5]=nd[2];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[3];
        a->pvccnode[count][2]=nd[4];
        a->pvccnode[count][3]=nd[5];
	
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[4];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[6];
        a->pvccnode[count][3]=nd[5];
		
        a->pvccnode[count][4]=nd[7];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[5];
        a->pvccnode[count][1]=nd[6];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[7];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		// cell 48
        if(pt[2]==0 && pt[3]==0 && cl[3]<=0 && cl[6]<=0 && cl[7]<=0)
        {
		// wedge 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[5];
        a->pvccnode[count][2]=nd[8];


        a->pvccnode[count][3]=nd[1];
        a->pvccnode[count][4]=nd[6];
        a->pvccnode[count][5]=nd[7];

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=6;

        ++count;
		
		// pyramid 1
        a->pvccnode[count][0]=nd[0];
        a->pvccnode[count][1]=nd[8];
        a->pvccnode[count][2]=nd[7];
        a->pvccnode[count][3]=nd[1];
	
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// pyramid 2
        a->pvccnode[count][0]=nd[2];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[8];
        a->pvccnode[count][3]=nd[3];
		
        a->pvccnode[count][4]=nd[4];
		
        a->pvccnode[count][5]=-1;
        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=5;

        ++count;
		
		// tetra 1
        a->pvccnode[count][0]=nd[1];
        a->pvccnode[count][1]=nd[7];
        a->pvccnode[count][2]=nd[2];


        a->pvccnode[count][3]=nd[4];
        a->pvccnode[count][4]=-1;
        a->pvccnode[count][5]=-1;

        a->pvccnode[count][6]=-1;
        a->pvccnode[count][7]=-1;

        a->ccedge[count]=4;

        ++count;
        }
		
		
		
        
}