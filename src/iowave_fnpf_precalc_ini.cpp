/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"


void iowave::fnpf_precalc_relax_ini(lexer *p, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    
    int dbcount=0;
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=0;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
        
    }
    
    // FI ------------------------------------------------
    FLOOP
    {
		dg = distgen(p); 
        db = distbeach(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ppt_count;

		}
        
        if(p->B99==1||p->B99==2)
		{
            if(db<1.0e20)
            ++dbcount;
        }
    }	

// ETA ------------------------------------------------
    SLICEBASELOOP
    {
		dg = distgen(p); 

		// Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++ept_count;

		}
    }	
    
    // precalc array alloc
    Fival.resize(ppt_count);
    Fifsfval.resize(ept_count);
    Fifsfval0.resize(ept_count);
    Fifsfval1.resize(ept_count);

    if(p->B89==1) 
    {
        etaval_S_sin.resize(ept_count,std::vector<double>(wave_comp));
        Fifsfval_S_sin.resize(ept_count,std::vector<double>(wave_comp));

        etaval_S_cos.resize(ept_count,std::vector<double>(wave_comp));
        Fifsfval_S_cos.resize(ept_count,std::vector<double>(wave_comp));

        etaval_T_sin.resize(wave_comp);
        Fifsfval_T_sin.resize(wave_comp);

        etaval_T_cos.resize(wave_comp);
        Fifsfval_T_cos.resize(wave_comp);
    }
}

void iowave::fnpf_precalc_dirichlet_ini(lexer *p, ghostcell *pgc)
{    
    int dbcount=0;
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=dbcount=0;
    
    
    // FI ------------------------------------------------
    FLOOP
    {
        db = distbeach(p); 
        
        if(p->B99==1||p->B99==2)
		{
            if(db<1.0e20)
            ++dbcount;
        }
    }	
    
    GCSLIN
    {
    i=p->gcslin[n].i;
    j=p->gcslin[n].j;
        FKLOOP
        ++ept_count;
    }
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
    }
      
    // precalc array alloc
    uval.resize(upt_count);
    Fival.resize(ppt_count);
    Uinval.resize(ppt_count);
    Fifsfval.resize(ept_count);
    Fifsfval0.resize(ept_count);
    Fifsfval1.resize(ept_count);

    if(p->B89==1) 
    {
        etaval_S_sin.resize(ept_count,std::vector<double>(wave_comp));
        Fival_S_sin.resize(ppt_count,std::vector<float>(wave_comp));
        Fifsfval_S_sin.resize(ept_count,std::vector<double>(wave_comp));
        uval_S_sin.resize(ppt_count,std::vector<float>(wave_comp));

        etaval_S_cos.resize(ept_count,std::vector<double>(wave_comp));
        Fival_S_cos.resize(ppt_count,std::vector<float>(wave_comp));
        Fifsfval_S_cos.resize(ept_count,std::vector<double>(wave_comp));
        uval_S_cos.resize(ppt_count,std::vector<float>(wave_comp));

        etaval_T_sin.resize(wave_comp);
        Fival_T_sin.resize(wave_comp);
        Fifsfval_T_sin.resize(wave_comp);
        uval_T_sin.resize(wave_comp);

        etaval_T_cos.resize(wave_comp);
        Fival_T_cos.resize(wave_comp);
        Fifsfval_T_cos.resize(wave_comp);
        uval_T_cos.resize(wave_comp);
    }
}
