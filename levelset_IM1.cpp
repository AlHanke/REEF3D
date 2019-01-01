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

#include"levelset_IM1.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reini.h"
#include"particlecorr.h"
#include"picard.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_fsf_entrain.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"

levelset_IM1::levelset_IM1(lexer* p, fdm *a, ghostcell* pgc, heat *&pheat, concentration *&pconc, field &ls):gradient(p), phin(p)
{

	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

	if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf(p,a,pgc);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_comp(p,a,pgc);
	
	if(p->F30>0 && p->H10>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);
	
	if(p->F30>0 && p->C10>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);
	
	if(p->F30>0 && p->F101>0 && p->F300==0 && p->W90==0)
	pupdate = new fluid_update_fsf_entrain(p,a,pgc,pconc);
	
	if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
	pupdate = new fluid_update_rheology(p,a,pgc);
	
	if(p->F300>0)
	pupdate = new fluid_update_void();


	if(p->F46==2)
	ppicard = new picard_f(p);

	if(p->F46==3)
	ppicard = new picard_lsm(p);

	if(p->F46!=2 && p->F46!=3)
	ppicard = new picard_void(p);
	
	ltimesave(p,a,ls);
}

levelset_IM1::~levelset_IM1()
{
}

void levelset_IM1::start(fdm* a,lexer* p, convection* pconvec,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particlecorr* ppart, field &ls)
{

    pflow->fsfinflow(p,a,pgc);
    ppicard->volcalc(p,a,pgc,ls);

    starttime=pgc->timer();
    clearrhs(p,a);
	pconvec->start(p,a,ls,4,a->u,a->v,a->w);
    timesource(p,a,pflow);
	psolv->start(p,a,pgc,ls,a->xvec,a->rhsvec,4,gcval_phi,p->F19);
	pflow->phi_relax(p,pgc,ls);
	pgc->start4(p,ls,gcval_phi);

    if(innercounter==p->N50-1 || p->N52==0)
    {
    ppart->start(p,a,pgc,pflow);
    pgc->start4(p,ls,gcval_phi);
    }
	
	p->lsmtime+=pgc->timer()-starttime;
	p->lsmiter+=p->solveriter;

    if(p->count%p->F41==0)
	preini->start(a,p,ls,pgc,pflow);

    ppicard->correct_ls(p,a,pgc,ls);
	ppart->picardmove(p,a,pgc);

	pflow->periodic(ls,p);
	pupdate->start(p,a,pgc);
    
	if(p->mpirank==0 && (innercounter==p->N50-1 || p->N52==0) && (p->count%p->P12==0))
	{

	cout<<"lsmiter: "<<p->lsmiter<<"  lsmtime: "<<setprecision(3)<<p->lsmtime<<endl;
	p->lsmiter=0;
	p->lsmtime=0.0;
	}
}

void levelset_IM1::timesource(lexer* p, fdm* a, ioflow *pflow)
{
    count=0;
    LOOP
    {
        a->M.p[count]+= 1.0/PDT;

        a->M.p[count]/= p->N56;

        a->rhsvec.V[count] += a->L(i,j,k) + phin(i,j,k)/PDT + a->M.p[count]*phin(i,j,k)*(1.0/p->N56-1.0);

	++count;
    }
}

void levelset_IM1::ltimesave(lexer* p, fdm *a, field &ls)
{
    LOOP
    phin(i,j,k)=ls(i,j,k);
}

void levelset_IM1::update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    pupdate->start(p,a,pgc);
}

void levelset_IM1::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	a->L(i,j,k)=0.0;
	++count;
    }
}

