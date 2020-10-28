/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include<iostream>
#include<fstream>
#include"6DOF_fsi.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void sixdof_fsi::print_ini(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank==0 && p->P14==1)
    {
        mkdir("./REEF3D_CFD_6DOF_STL", 0777);
        mkdir("./REEF3D_CFD_6DOF", 0777);
    }
	
    ofstream print;
    
	if(p->P14==0)
	print.open("REEF3D_6DOF_position.dat");
	if(p->P14==1)
	print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_position.dat");
	
	print<<"time \t XG \t YG \t ZG \t Phi \t Theta \t Psi"<<endl;
	
	print.close();
    
	
	if(p->P14==0)
	print.open("REEF3D_6DOF_velocity.dat");
	if(p->P14==1)
	print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_velocity.dat");
	
	print<<"time \t Ue [m/s] \t Ve [m/s] \t We [m/s] \t Pe [rad/s] \t Qe [rad/s] \t Re [rad/s]"<<endl;
	
    print.close();
    

	if(p->P14==0)
	print.open("REEF3D_6DOF_forces.dat");
	if(p->P14==1)
	print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_forces.dat");
	
	print<<"time \t Fx \t Fy \t Fz \t Mx \t My \t Mz \t cd \t cq \t cl"<<endl;
	
    print.close();    
	

    if(p->P14==0)
	print.open("REEF3D_6DOF_surface_forces.dat");
	if(p->P14==1)
	print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_surface_forces.dat");
	
	print<<"time \t Fx \t Fy \t Fz \t Mx \t My \t Mz \t cd \t cq \t cl"<<endl;
	
    print.close();    
}


void sixdof_fsi::print_stl(lexer *p, fdm *a, ghostcell *pgc)
{
	int num=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;

    if(p->mpirank==0 && (((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || p->count==0))
    {
        printtime+=p->P30;
        
        char path[100];
        
        if(p->P14==1)
        {
            if(num<10)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-00000%d.stl",num);

            if(num<100&&num>9)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-0000%d.stl",num);

            if(num<1000&&num>99)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-000%d.stl",num);

            if(num<10000&&num>999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-00%d.stl",num);

            if(num<100000&&num>9999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-0%d.stl",num);

            if(num>99999)
            sprintf(path,"./REEF3D_CFD_6DOF_STL/REEF3D-6DOF-%d.stl",num);
        }

        ofstream result;
        result.open(path, ios::binary);

        
        result<<"solid"<<" "<<"ascii"<<endl;
        
        double zero=0.0;
        
        for(n=0; n<tricount; ++n)
        {
        result<<" facet normal "<<zero<<" "<<zero<<" "<<zero<<endl;
        result<<"  outer loop"<<endl;
        result<<"   vertex "<<tri_x[n][0]<<" "<<tri_y[n][0]<<" "<<tri_z[n][0]<<endl;
        result<<"   vertex "<<tri_x[n][1]<<" "<<tri_y[n][1]<<" "<<tri_z[n][1]<<endl;
        result<<"   vertex "<<tri_x[n][2]<<" "<<tri_y[n][2]<<" "<<tri_z[n][2]<<endl;
        result<<"  endloop"<<endl;
        result<<" endfacet"<<endl;
        }
        
        result<<"endsolid"<<endl;
        

        result.close();

        ++p->printcount_sixdof;	
    }
}


void sixdof_fsi::print_parameter(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->mpirank == 0 && p->count%p->X19==0)
    {
        ofstream print;
        
        if(p->P14==0)
        print.open("REEF3D_6DOF_position.dat", ofstream::app);
        if(p->P14==1)
        print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_position.dat", ofstream::app);
        
        //print<<p->simtime<<" \t "<<c_(0)<<" \t "<<c_(1)<<" \t "<<c_(2)<<" \t "<<phi*(180.0/PI)<<" \t "<<theta*(180.0/PI)<<" \t "<<psi*(180.0/PI)<<endl;
        print<<p->simtime<<" \t "<<p->xg<<" \t "<<p->yg<<" \t "<<p->zg<<" \t "<<phi*(180/PI)<<" \t "<<theta*(180/PI)<<" \t "<<psi*(180/PI)<<endl;

        print.close();
        
        
        if(p->P14==0)
        print.open("REEF3D_6DOF_velocity.dat", ofstream::app);
        if(p->P14==1)
        print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_velocity.dat", ofstream::app);
        
        print<<p->simtime<<" \t "<<p->ufbi<<" \t "<<p->vfbi<<" \t "<<p->wfbi<<" \t "<<p->pfbi<<" \t "<<p->qfbi<<" \t "<<p->rfbi<<endl;
        
        print.close();
/*
        if(p->P14==0)
        print.open("REEF3D_6DOF_forces.dat", ofstream::app);
        if(p->P14==1)
        print.open("./REEF3D_CFD_6DOF/REEF3D_6DOF_forces.dat", ofstream::app);

        print<<p->simtime<<" \t "<<Ffb_(0)<<" \t "<<Ffb_(1)<<" \t "<<Ffb_(2)<<" \t "<<Mfb_(0)<<" \t "<<Mfb_(1)<<" \t "<<Mfb_(2)<<" \t "<<0.0<<" \t "<<0.0<<" \t "<<0.0<<endl;  
        print.close();
*/
    }
}
