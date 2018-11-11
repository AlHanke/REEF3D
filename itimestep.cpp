/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"itimestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

itimestep::itimestep(lexer* p):epsi(1.0e-19),maxtimestep(p->N49),c0_orig(p->N47)
{
}

itimestep::~itimestep()
{
}

void itimestep::start(fdm *a, lexer *p,ghostcell *pgc, turbulence *pturb)
{
    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->dt_old=p->dt;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
	sqd=1.0/(p->DXM*p->DXM);

// maximum velocities

	ULOOP
	p->umax=MAX(p->umax,fabs(a->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	VLOOP
	p->vmax=MAX(p->vmax,fabs(a->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	WLOOP
	p->wmax=MAX(p->wmax,fabs(a->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);
	

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<endl;
    }
	
	p->umax=MAX(p->umax,p->ufbmax);
	p->vmax=MAX(p->vmax,p->vfbmax);
	p->wmax=MAX(p->wmax,p->wfbmax);

// maximum viscosity
	LOOP
	p->viscmax=MAX(p->viscmax, a->visc(i,j,k)+a->eddyv(i,j,k));

	p->viscmax=pgc->globalmax(p->viscmax);

    if(p->mpirank==0)
	cout<<"viscmax: "<<p->viscmax<<endl;
	//----kin
	LOOP
	p->kinmax=MAX(p->kinmax,pturb->kinval(i,j,k));

	p->kinmax=pgc->globalmax(p->kinmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kinmax: "<<p->kinmax<<endl;

	//---eps
    LOOP
	p->epsmax=MAX(p->epsmax,pturb->epsval(i,j,k));

	p->epsmax=pgc->globalmax(p->epsmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"epsmax: "<<p->epsmax<<endl;

	//---press
    LOOP
    {
	p->pressmax=MAX(p->pressmax,a->press(i,j,k));
	p->pressmin=MIN(p->pressmin,a->press(i,j,k));
    }

	p->pressmax=pgc->globalmax(p->pressmax);
	p->pressmin=pgc->globalmin(p->pressmin);


//  increasing N4/
    if(p->N57_1>0 && p->N57_2>0)
    {
        if(p->count>=p->N57_1 && p->count<p->N57_2)
        p->N47=c0_orig + (p->N57_3-c0_orig)/(double(p->N57_2)-double(p->N57_1))*(p->count-p->N57_1);

        if(p->count>=p->N57_2)
        p->N47=p->N57_3;
    }

//  connect residuals to time marching
    if(p->N58==1)
    {
    p->N53/=(1.0+p->N47);
	p->N54/=(1.0+p->N47);
	p->N55/=(1.0+p->N47);
	p->N56/=(1.0+p->N47);
    }

// maximum reynolds stress source term
	visccrit=p->viscmax*(6.0/pow(p->DXM,2.0));
    

    cu=cv=cw=1.0e10;
    ULOOP
    {
    dx = MIN3(p->DXP[IP],p->DYN[JP],p->DZN[KP]);

	cu = MIN(cu, 2.0/((fabs(p->umax)/dx)));
    }

    VLOOP
    {
    dx = MIN3(p->DXN[IP],p->DYP[JP],p->DZN[KP]);
    
	cv = MIN(cv, 2.0/((fabs(p->vmax)/dx)));
    }
    
    WLOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZP[KP]);

	cw = MIN(cw, 2.0/((fabs(p->wmax)/dx)));
    }
    

   	p->dt=p->N47*min(cu,cv,cw);
	p->dt=pgc->timesync(p->dt);

	p->dt=MIN(p->dt,maxtimestep);
	p->turbtimestep=p->dt;
	

	a->maxF=0.0;
	a->maxG=0.0;
	a->maxH=0.0;
}

void itimestep::ini(fdm* a, lexer* p,ghostcell* pgc)
{

	p->umax=p->vmax=p->wmax=p->viscmax=-1e19;
	p->umax=p->W10;
	
	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);
    
    p->umax=MAX(p->umax,2.0);

	p->dt=p->DXM/(p->umax+epsi);

	LOOP
	p->viscmax=MAX(p->viscmax, a->visc(i,j,k)+a->eddyv(i,j,k));

	visccrit=(p->viscmax*(6.0/pow(p->DXM,2.0)));

	cu=2.0/((p->umax/p->DXM+visccrit)+sqrt(pow(p->umax/p->DXM+visccrit,2.0)+(4.0*sqrt(fabs(a->gi) + fabs(a->gj) +fabs(a->gk)))/p->DXM));// + (8.0*p->maxkappa*p->W5)/(2.0*p->DXM*p->DXM*(p->W1+p->W3)));

    //  connect residuals to time marching
    if(p->N58==1)
    {
    p->N53/=(1.0+p->N47);
	p->N54/=(1.0+p->N47);
	p->N55/=(1.0+p->N47);
	p->N56/=(1.0+p->N47);
    }

	p->dt=p->N47*cu;
	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;

	p->maxkappa=0.0;
}

double itimestep::min(double val1,double val2,double val3)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini>val3)
	mini=val3;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double itimestep::min(double val1,double val2)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double itimestep::max(double val1,double val2,double val3)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<val3)
	maxi=val3;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

double itimestep::max(double val1,double val2)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

