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

#include"wave_lib_piston_eta.h"
#include"lexer.h"
#include<fstream>

wave_lib_piston_eta::wave_lib_piston_eta(lexer *p) : wave_lib_parameters(p)
{ 
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: piston_eta wavemaker theory";
    }
	
	timecount=0;
	
	read(p);
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));	
}

double wave_lib_piston_eta::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_piston_eta::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_piston_eta::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;

    if(p->wavetime<ts || p->wavetime>te || timecount>=ptnum-1)
	return 0.0;
	
	if(p->wavetime>eta[timecount+1][0])
	++timecount;
	
	vel = sqrt(9.81/wdt) * wave_eta(p,x,y);
    
    if(p->B110==1)
    {
    z+=p->wd;
    
    if(z<p->B110_zs || z>p->B110_ze)
    vel=0.0;
    }

    return vel;
}

double wave_lib_piston_eta::wave_w(lexer *, double, double, double)
{
    return 0.0;
}

double wave_lib_piston_eta::wave_eta(lexer *p, double, double)
{
    double val=0.0;
    
    if(p->wavetime<ts || p->wavetime>te || timecount>=ptnum-1)
	return 0.0;
    
    val =  ((eta[timecount+1][1]-eta[timecount][1])/(eta[timecount+1][0]-eta[timecount][0]))
            *((p->wavetime)-eta[timecount][0]) + eta[timecount][1];
	
    return val;
}

double wave_lib_piston_eta::wave_fi(lexer *, double, double, double)
{
    return 0.0;
}

void wave_lib_piston_eta::parameters(lexer*)
{

}

void wave_lib_piston_eta::read(lexer* p)
{
	double val0,val1;
	int count;
	
	char name[] = "wavemaker_eta.dat";

// open file------------
	ifstream file(name, ios_base::in);
	
	if(!file)
	{
		cout<<endl<<("no 'wavemaker_eta.dat' file found")<<endl<<endl;

	}
	
	count=0;
	while(!file.eof())
	{
	file>>val0>>val1;
	if(val0>=p->B117)
	++count;
	}
	
	file.close();

    
    ptnum=count;
	
	p->Darray(eta,ptnum,2);
	
	file.open (name, ios_base::in);
	
	count=0;
	while(!file.eof())
	{
	
	file>>val0>>val1;
	
	if(val0>=p->B117)
	{
	eta[count][0] = val0-p->B117;
	eta[count][1] = val1;
	++count;
	}
	}
	
	ts = eta[0][0];
	te = eta[ptnum-1][0];
	
}

void wave_lib_piston_eta::wave_prestep(lexer*)
{
}
