/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/


#ifndef SEDIMENT_FDM_H_
#define SEDIMENT_FDM_H_

#include"sliceint4.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"field4a.h"

using namespace std;

class sediment_fdm
{
public:
    sediment_fdm(lexer*);
    virtual ~sediment_fdm();
    
    slice1 P; ///< u vel on bed
    slice2 Q; ///< v vel on bed
    
    slice4 bedzh; ///< bed elevation
    slice4 bedzh0; ///< bed elevation, initial
    slice4 bedch; ///< bed elevation, change
    slice4 bedsole; ///< unused
    slice4 vz; ///< bed change vel
    slice4 dh; ///< bed change vel
    slice4 reduce; ///< critical hear stress reduction factor
    slice4 ks; ///< bed roughnes
    
    slice4 tau_eff; ///< effective shear stress
    slice4 tau_crit; ///< critical shear stress
    slice4 shearvel_eff; ///< effective shear velocity
    slice4 shearvel_crit; ///< critical shear velocity
    slice4 shields_eff; ///< effective shields number
    slice4 shields_crit; ///< critical shields number
    
    slice4 alpha; ///< bed angle
    slice4 teta; ///< bed angle
    slice4 gamma; ///< bed angle
    slice4 beta; ///< bed angle
    slice4 phi; ///< bed angle
    sliceint4 active; ///< flag for active bed
    
    sliceint4 bedk; ///< k index of cell containing bed
    slice4 slideflag; ///< flag for bed slide
    
    slice4 qb; ///< bed load - exner
    slice4 qbe; ///< bed load - bedload forumla
    slice4 cbe; ///< bedconc ?
    slice4 cb; ///< bedconc ?
    slice4 cbn; ///< bedconc ?
    slice4 conc; ///< bedconc ?
    
    slice4 waterlevel;  ///< waterdepth - unused
    slice4 guard; ///< ?
    
    double ws; ///< single particle sedimenation velocity
};;

#endif