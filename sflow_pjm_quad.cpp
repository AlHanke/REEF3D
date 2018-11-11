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

#include"sflow_pjm_quad.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver2D.h"
#include"momentum.h"
#include"ioflow.h"
#include"sflow_weno_hj.h"
#include"sflow_gradient_weno.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HXIMJ (fabs(b->hx(i-1,j))>1.0e-20?b->hx(i-1,j):1.0e20)

#define HXP (0.5*(HXIJ + HXIMJ))
#define HYP (0.5*(HYIJ + HYIJM))

#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
#define HYIJM (fabs(b->hy(i,j-1))>1.0e-20?b->hy(i,j-1):1.0e20)

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

#define HPIP (fabs(b->hp(i+1,j))>1.0e-20?b->hp(i+1,j):1.0e20)
#define HPJP (fabs(b->hp(i,j+1))>1.0e-20?b->hp(i,j+1):1.0e20)

#define HPIM (fabs(b->hp(i-1,j))>1.0e-20?b->hp(i-1,j):1.0e20)
#define HPJM (fabs(b->hp(i,j-1))>1.0e-20?b->hp(i,j-1):1.0e20)

#define HPXP (0.5*(HP + HPIP))
#define HPYP (0.5*(HP + HPJP))

#define HPXM (0.5*(HP + HPIM))
#define HPYM (0.5*(HP + HPJM))
 
sflow_pjm_quad::sflow_pjm_quad(lexer* p, fdm2D *b) : phi4(p)
{
    if(p->B76==0)
    gcval_press=40;  

    if(p->B76==1)
    gcval_press=41;

    if(p->B76==2)
    gcval_press=42;

    if(p->B76==3)
    gcval_press=43;
	
	if(p->B76==4) 
    gcval_press=44;
	
	if(p->B76==5) 
    gcval_press=45;
	
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
	disc = new sflow_weno_hj(p);
    
    pgrad = new sflow_gradient_weno(p);
    
    
    
    wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->dx;
	
}

sflow_pjm_quad::~sflow_pjm_quad()
{
}

void sflow_pjm_quad::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, slice &P, slice &Q, slice &Pn, slice &Qn, slice &ws, slice &eta, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
	
	starttime=pgc->timer();
	
    quad_calc(p,b,P,Q,Pn,Qn,alpha);
    
    rhs(p,b,P,Q,ws,alpha);
	pgc->gcsl_start4(p,b->press,gcval_press);
	
    pgc->gcsl_start4(p,phi4,1);
    
    poisson(p,b,alpha);
	
        solvtime=pgc->timer();

    psolv->start(p,b,pgc,b->press,b->xvec,b->rhsvec,4,gcval_press,p->N44);
	
        p->poissontime=pgc->timer()-solvtime;
  
    pflow->pm_relax(p,pgc,b->press);
	pgc->gcsl_start4(p,b->press,gcval_press);
	b->press.ggcpol(p);
    
    
	ucorr(p,b,P,eta,alpha);
	vcorr(p,b,Q,eta,alpha);
    wcorr(p,b, alpha, P, Q, ws);
	pgc->gcsl_start4(p,ws,12);

    p->poissoniter=p->solveriter;

	ptime=pgc->timer()-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0) && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  solvtime: "<<setprecision(3)<<p->poissontime<<"  ptime: "<<setprecision(3)<<ptime<<endl;
}

void sflow_pjm_quad::ucorr(lexer* p, fdm2D* b, slice& P, slice &eta, double alpha)
{	
	SLICELOOP1
    if(b->breaking(i,j)==0 && b->breaking(i+1,j)==0)
	P(i,j) -= alpha*p->dt*(((b->press(i+1,j)-b->press(i,j))/(p->dx*p->W1)));

    if(p->A222==1)           
    SLICELOOP1
    if(b->breaking(i,j)==0 && b->breaking(i+1,j)==0)
	P(i,j) += alpha*p->dt*(0.75*(b->press(i+1,j)+b->press(i,j))*((b->depth(i+1,j)-b->depth(i,j))
                            /(p->dx*HPXP*p->W1))
    
                        + 0.125*(phi4(i+1,j)+phi4(i,j))*((b->depth(i+1,j)-b->depth(i,j))/(p->dx)));
}

void sflow_pjm_quad::vcorr(lexer* p, fdm2D* b, slice& Q, slice &eta, double alpha)
{	
	SLICELOOP2
    if(b->breaking(i,j)==0 && b->breaking(i,j+1)==0)
	Q(i,j) -= alpha*p->dt*(((b->press(i,j+1)-b->press(i,j))/(p->dx*p->W1)));
                
    if(p->A222==1)
    SLICELOOP2
    if(b->breaking(i,j)==0 && b->breaking(i,j+1)==0)
	Q(i,j) += alpha*p->dt*(0.75*(b->press(i,j+1)+b->press(i,j))*((b->depth(i,j+1)-b->depth(i,j))
                            /(p->dx*HPYP*p->W1))
                            
                        + 0.125*(phi4(i,j+1)+phi4(i,j))*((b->depth(i,j+1)-b->depth(i,j))/(p->dx)));
}

void sflow_pjm_quad::wcorr(lexer* p, fdm2D* b, double alpha, slice &P, slice &Q, slice &ws)
{	    
    SLICELOOP4
    if(b->breaking(i,j)==0)
	 ws(i,j) += p->dt*alpha*(1.5*b->press(i,j)/(HP*p->W1)  +  0.25*phi4(i,j));
}

void sflow_pjm_quad::wcalc(lexer* p, fdm2D* b,double alpha, slice &P, slice &Q, slice &ws)
{	
}

void sflow_pjm_quad::rhs(lexer *p, fdm2D* b, slice &P, slice &Q, slice &ws, double alpha)
{
    NSLICELOOP4
	b->rhsvec.V[n]=0.0;

    count=0;
    SLICELOOP4
    {
    b->rhsvec.V[count] =    - ((P(i,j) - P(i-1,j))*b->hp(i,j)
                             + (Q(i,j) - Q(i,j-1))*b->hp(i,j))/(alpha*p->dt*p->dx)
                           
                            - 2.0*
							(
								ws(i,j) 
								+ 0.25*(P(i,j)+P(i-1,j))*(b->depth(i+1,j)-b->depth(i-1,j))/p->dx
								+ 0.25*(Q(i,j)+Q(i,j-1))*(b->depth(i,j+1)-b->depth(i,j-1))/p->dx
                           )
							/(alpha*p->dt);
                           
    ++count;
    
    
    
    }
    
    if(p->A221<2)
    SLICELOOP4
    b->press(i,j)=0.0;
}

void sflow_pjm_quad::poisson(lexer*p, fdm2D* b, double alpha)
{
    sqd = (1.0/(p->dx*p->dx*p->W1));
        
     n=0;
    SLICELOOP4
	{
	b->M.p[n]  =  (b->hp(i,j)*sqd + b->hp(i,j)*sqd)*p->x_dir
               +  (b->hp(i,j)*sqd + b->hp(i,j)*sqd)*p->y_dir + 2.0/(HP*p->W1);     

               

   	b->M.n[n] = -b->hp(i,j)*sqd*p->x_dir;
	b->M.s[n] = -b->hp(i,j)*sqd*p->x_dir;

	b->M.w[n] = -b->hp(i,j)*sqd*p->y_dir;
	b->M.e[n] = -b->hp(i,j)*sqd*p->y_dir;

	++n;
	}
	
	
    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0)
		{
		b->rhsvec.V[n] -= b->M.s[n]*b->press(i-1,j);
		b->M.s[n] = 0.0;
		}
		
		if(p->flagslice4[Ip1J]<0)
		{
		b->rhsvec.V[n] -= b->M.n[n]*b->press(i+1,j);
		b->M.n[n] = 0.0;
		}
		
		if(p->flagslice4[IJm1]<0)
		{
		b->rhsvec.V[n] -= b->M.e[n]*b->press(i,j-1);
		b->M.e[n] = 0.0;
		}
		
		if(p->flagslice4[IJp1]<0)
		{
		b->rhsvec.V[n] -= b->M.w[n]*b->press(i,j+1);
		b->M.w[n] = 0.0;
		}
		
	++n;
	}
    
    n=0;
    SLICELOOP4
	{
        if(b->wet4(i,j)==0 || b->breaking(i,j)==1)
        {
        b->M.p[n]  = 1.0;

        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;

        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        
        b->rhsvec.V[n]=0.0;
        }
	++n;
	}
}

void sflow_pjm_quad::upgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
        if(p->A221==1)
        SLICELOOP1
        b->F(i,j) -= fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) 
                                 - p->A223*eta(i,j) -  (1.0-p->A223)*eta_n(i,j) )/(p->dx); 
                                
        if(p->A221==2)           
        SLICELOOP1
        if(b->breaking(i,j)==0 && b->breaking(i+1,j)==0)
        b->F(i,j) += -(((b->press(i+1,j)-b->press(i,j))/(p->dx*p->W1)));
}

void sflow_pjm_quad::vpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
        if(p->A221==1)
        SLICELOOP2
        b->G(i,j) -= fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) 
                                 - p->A223*eta(i,j) -  (1.0-p->A223)*eta_n(i,j) )/(p->dx); 
        
        if(p->A221==2)
        SLICELOOP2
        if(b->breaking(i,j)==0 && b->breaking(i,j+1)==0)
        b->G(i,j) +=   - (((b->press(i,j+1)-b->press(i,j))/(p->dx*p->W1)));
}

void sflow_pjm_quad::quad_calc(lexer *p,fdm2D *b,slice &P, slice &Q, slice &Pn, slice &Qn, double alpha)
{
    /*
    SLICELOOP1
    {
    phi1(i,j) = - ((b->depth(i+1,j)-b->depth(i,j))/p->dx) 
    
                                * ((P(i,j)-Pn(i,j))/(alpha*p->dt)
    
                                + MAX(0.0,P(i,j))*((P(i,j)-P(i-1,j))/p->dx) + MIN(0.0,P(i,j))*((P(i+1,j)-P(i,j))/p->dx)
                                + MAX(0.0,0.5*(Q(i,j-1)+Q(i+1,j-1)))*((P(i,j)-P(i,j-1))/p->dx) + MIN(0.0,0.5*(Q(i,j)+Q(i+1,j)))*((P(i,j+1)-P(i,j))/p->dx))
                                                          
                - pow(b->P(i,j),2.0)*((b->depth(i+2,j) - b->depth(i+1,j) - b->depth(i,j) + b->depth(i-1,j))/(2.0*p->dx*p->dx));
    
    }
    
    SLICELOOP2
    {
    phi2(i,j) = - ((b->depth(i,j+1)-b->depth(i,j))/p->dx) 
                                * ((Q(i,j)-Qn(i,j))/(alpha*p->dt)
    
                                + MAX(0.0,0.5*(P(i-1,j)+P(i-1,j+1)))*((Q(i,j)-Q(i-1,j))/p->dx) + MIN(0.0,0.5*(P(i,j)+P(i,j+1)))*((Q(i+1,j)-Q(i,j))/p->dx)
                                + MAX(0.0,Q(i,j))*((Q(i,j)-Q(i,j-1))/p->dx) + MIN(0.0,Q(i,j))*((Q(i,j+1)-Q(i,j))/p->dx))
                                                          
                - pow(b->Q(i,j),2.0)*((b->depth(i,j+2) - b->depth(i,j+1) - b->depth(i,j) + b->depth(i,j-1))/(2.0*p->dx*p->dx));
    
    }*/
    
    double Pval,Pnval;
    double Qval,Qnval;
    
    
    SLICELOOP4
    {   
        Pval = 0.5*(P(i-1,j)+P(i,j));
        Pnval = 0.5*(Pn(i-1,j)+Pn(i,j));
        
        Qval = 0.5*(Q(i,j-1)+Q(i,j));
        Qnval = 0.5*(Qn(i,j-1)+Qn(i,j));
    

    phi4(i,j) = -(MAX(0.0,Pval)*((b->depth(i,j)-b->depth(i-1,j))/(p->dx))
                    
                + MIN(0.0,Pval)*((b->depth(i+1,j)-b->depth(i,j))/(p->dx)))

    
                                    * ((Pval-Pnval)/(alpha*p->dt)
                                                                    
                                    + Pval*((P(i,j)-P(i-1,j))/(p->dx)) 
                                    
                                    //+ Pval*(P(i-2,j) - 27.0*P(i-1,j) + 27.0*P(i,j) - P(i+1,j))/(24.0*p->dx) 
                                                         
                                    + Qval*((0.5*(P(i,j+1)+P(i-1,j+1))-0.5*(P(i,j-1)+P(i-1,j-1)))/(2.0*p->dx)))
                                                

                -(MAX(0.0,Qval)*((b->depth(i,j)-b->depth(i,j-1))/(p->dx))
                    
                + MIN(0.0,Qval)*((b->depth(i,j+1)-b->depth(i,j))/(p->dx)))
                
                                    * ((Qval-Qnval)/(alpha*p->dt)    
     
                                    + Pval*((0.5*(Q(i+1,j)+Q(i+1,j-1))-0.5*(Q(i-1,j)+Q(i-1,j-1)))/(2.0*p->dx))
                                                         
                                    + Qval*((Q(i,j-1)-Q(i,j))/p->dx))
                                            
                                                          
                - pow(Pval,2.0)*((b->depth(i+1,j) - 2.0*b->depth(i,j) + b->depth(i-1,j))/(p->dx*p->dx));
                
                - pow(Qval,2.0)*((b->depth(i,j+1) - 2.0*b->depth(i,j) + b->depth(i,j-1))/(p->dx*p->dx));
       
    }
}

void sflow_pjm_quad::wpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{	    
}
