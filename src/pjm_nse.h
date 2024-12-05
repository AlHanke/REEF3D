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

#ifndef PJM_NSE_H_
#define PJM_NSE_H_

#include"pressure.h"
#include"increment.h"

class heat;
class concentration;
class density;

class pjm_nse : public pressure, public increment
{

public:

	pjm_nse(lexer* , fdm*, heat*&, concentration*&);
	virtual ~pjm_nse() = default;

	virtual void start(lexer*,fdm*,ghostcell*,poisson*, solver*, ioflow*, field&, field&, field&,double);
    virtual void ini(lexer*,fdm*,ghostcell*);
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void ucorr(lexer*p,fdm*,field&,double);
	virtual void vcorr(lexer*p,fdm*,field&,double);
	virtual void wcorr(lexer*p,fdm*,field&,double);
	virtual void upgrad(lexer*,fdm*,slice&,slice&);
	virtual void vpgrad(lexer*,fdm*,slice&,slice&);
    virtual void wpgrad(lexer*,fdm*,slice&,slice&);

private:
	double starttime,endtime;
    double teta;
	int count, gcval_press, check;
	int gcval_u, gcval_v, gcval_w;
	
	void debug(lexer*,fdm*);
    
    density *pd;
};



#endif

