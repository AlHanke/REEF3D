/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include"mooring_Catenary.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void mooring_Catenary::print(lexer *p)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;
	
	if
	(
		p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  
		|| (p->simtime>printtime && p->P30>0.0)   
		|| p->count==0)
	)
	{
		printtime+=p->P30;
		
		if(p->P14==1)
		{
			if(num<10)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-00000%d.vtk",line,num);

			if(num<100&&num>9)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-0000%d.vtk",line,num);

			if(num<1000&&num>99)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-000%d.vtk",line,num);

			if(num<10000&&num>999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-00%d.vtk",line,num);

			if(num<100000&&num>9999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-0%d.vtk",line,num);

			if(num>99999)
			sprintf(name,"./REEF3D_6DOF_Mooring/REEF3D-Mooring-%d-%d.vtk",line,num);
		}
		
		// Reconstruct line
		buildLine(p);
		

		// Print results
		ofstream result;
		result.open(name, ios::binary);
		
		result << "# vtk DataFile Version 2.0" << endl;
		result << "Mooring line " << line << endl;
		result << "ASCII \nDATASET UNSTRUCTURED_GRID" << endl;
		result << "POINTS " << H << " float" <<endl;

		for (int n = 0; n < H; ++n)
		{
			result<<x[n]<<" "<<y[n]<<" "<<z[n]<<endl;
		}
		
		result << "\nCELLS " << H-1 << " " << (H-1)*3 <<endl;	
		
		for(int n = 0; n < (H-1); ++n)
		{
			result<<"2 "<< n << " " << n+1 << endl;
		}
		
		result << "\nCELL_TYPES " << H-1 << endl;	
		
		for(int n = 0; n < (H-1); ++n)
		{
			result<<"3"<<endl;
		}	

		result<<"\nPOINT_DATA " << H <<endl;
		result<<"SCALARS Tension float 1 \nLOOKUP_TABLE default"<<endl;
		
		for(int n = 0; n < H; ++n)
		{
			result<<T[n]<<endl;
		}
		
		result.close();


		eTout<<p->simtime<<" \t "<<T[H-1]<<endl;	
	}
}


void mooring_Catenary::buildLine(lexer *p)
{
	double d_xy,dl, segLen, alpha;

	if (EA == 0.0)
	{
		dl = L - lms;
		segLen = (dxy - dl)/(H-1);			
		alpha = atan(dy/dx);

		x[0] = xs;
		y[0] = ys;
		z[0] = zs;
		T[0] = fabs(FH);
			
		z[1] = zs;
		T[1] = fabs(FH);
			
		if (dx > 0)
		{
			x[1] = xs + dl*cos(alpha);
		}
		else
		{
			x[1] = xs - dl*cos(alpha);
		}
		if (dy > 0)
		{
			y[1] = ys + dl*sin(alpha);
		}
		else
		{
			y[1] = ys - dl*sin(alpha);
		}
			
		for (int cnt = 2; cnt < H; cnt++)
		{
			if (dx > 0)
			{
				x[cnt] = x[1] + cnt*segLen*cos(alpha);
			}
			else
			{
				x[cnt] = x[1] - cnt*segLen*cos(alpha);
			}
			if (dy > 0)
			{
				y[cnt] = y[1] + cnt*segLen*sin(alpha);
			}
			else
			{
				y[cnt] = y[1] - cnt*segLen*sin(alpha);
			}		
			
			z[cnt] = fabs(FH)/w*(cosh(w/fabs(FH)*fabs(cnt*segLen)) - 1.0);
			T[cnt] = fabs(FH) + w*z[cnt];
		}	
	}
	else
	{
		segLen = L/(H-1);
		alpha = atan(dy/dx);

		for (int cnt = 0; cnt < H; cnt++)
		{
			if (segLen*cnt <= lms)
			{
				z[cnt] = zs;
				T[cnt] = fabs(FH);
						
				d_xy = segLen*cnt*(1.0 + FH/EA);
						
				if (dx > 0)
				{
					x[cnt] = xs + d_xy*cos(alpha);
				}
				else
				{
					x[cnt] = xs - d_xy*cos(alpha);
				}
					
				if (dy > 0)
				{
					y[cnt] = ys + d_xy*sin(alpha);
				}
				else
				{
					y[cnt] = ys - d_xy*sin(alpha);
				}
			}
			else
			{
				z[cnt] = zs + FH/w*(sqrt(1+(w*(segLen*cnt-lms)/FH)*(w*(segLen*cnt-lms)/FH)) - 1) + w*(segLen*cnt-lms)*(segLen*cnt-lms)/(2*EA);		
				T[cnt] = sqrt(FH*FH + (w*(segLen*cnt - lms))*(w*(segLen*cnt - lms)));
						
				d_xy = FH/w*log(w*(segLen*cnt-lms)/FH + sqrt(1+(w*(segLen*cnt-lms)/FH)*(w*(segLen*cnt-lms)/FH))) + FH*segLen*cnt/EA;
					
				if (dx > 0)
				{
					x[cnt] = xs + lms*cos(alpha) + d_xy*cos(alpha);
				}
				else
				{
					x[cnt] = xs - lms*cos(alpha) - d_xy*cos(alpha);
				}
						
				if (dy > 0)
				{
					y[cnt] = ys + lms*sin(alpha) + d_xy*sin(alpha);
				}
				else
				{
					y[cnt] = ys - lms*sin(alpha) - d_xy*sin(alpha);
				}					
			}
		}
	}
}