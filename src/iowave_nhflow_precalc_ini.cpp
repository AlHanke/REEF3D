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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::nhflow_precalc_relax_ini(lexer *p,fdm_nhf *d, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    
    upt_count=ept_count=0;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
    }
    
    // U ------------------------------------------------
    BASELOOP
    {
        dg = distgen(p);
        
        // Wave Generation
		if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            ++upt_count;
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
    
    //cout<<p->mpirank<<"  ept_count: "<<ept_count<<"  upt_count: "<<upt_count<<endl;
    
    // precalc array alloc
    uval.resize(upt_count);
    vval.resize(upt_count);
    wval.resize(upt_count);
    UHval.resize(upt_count);
    VHval.resize(upt_count);
    WHval.resize(upt_count);

    if(p->B89==1) 
    {
        uval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        vval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        wval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        etaval_S_sin.resize(ept_count,std::vector<double>(wave_comp));

        uval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        vval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        wval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        etaval_S_cos.resize(ept_count,std::vector<double>(wave_comp));

        uval_T_sin.resize(wave_comp);
        vval_T_sin.resize(wave_comp);
        wval_T_sin.resize(wave_comp);
        etaval_T_sin.resize(wave_comp);

        uval_T_cos.resize(wave_comp);
        vval_T_cos.resize(wave_comp);
        wval_T_cos.resize(wave_comp);
        etaval_T_cos.resize(wave_comp);
    }
}

void iowave::nhflow_precalc_dirichlet_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // count number of relax points
    // allocate double* array
    
    upt_count=ept_count = p->gcin_count;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
    }
  
    // precalc array alloc
    uval.resize(upt_count);
    vval.resize(upt_count);
    wval.resize(upt_count);
    UHval.resize(upt_count);
    VHval.resize(upt_count);
    WHval.resize(upt_count);

    if(p->B89==1) 
    {
        uval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        vval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        wval_S_sin.resize(upt_count,std::vector<float>(wave_comp));
        etaval_S_sin.resize(ept_count,std::vector<double>(wave_comp));

        uval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        vval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        wval_S_cos.resize(upt_count,std::vector<float>(wave_comp));
        etaval_S_cos.resize(ept_count,std::vector<double>(wave_comp));

        uval_T_sin.resize(wave_comp);
        vval_T_sin.resize(wave_comp);
        wval_T_sin.resize(wave_comp);
        etaval_T_sin.resize(wave_comp);

        uval_T_cos.resize(wave_comp);
        vval_T_cos.resize(wave_comp);
        wval_T_cos.resize(wave_comp);
        etaval_T_cos.resize(wave_comp);
    }
}
