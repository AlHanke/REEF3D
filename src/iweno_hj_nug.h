/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef IWENO_HJ_NUG_H_
#define IWENO_HJ_NUG_H_

#include"convection.h"
#include"weno_nug_func.h"

class flux;

class iweno_hj_nug : public convection, public weno_nug_func
{
public:
    iweno_hj_nug (lexer*);
    virtual ~iweno_hj_nug() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) final;

private:
    template<typename GenericField>
    void wenoloop1(lexer*,fdm*,const GenericField&,int,const GenericField&,const GenericField&,const GenericField&);
    template<typename GenericField>
    void wenoloop2(lexer*,fdm*,const GenericField&,int,const GenericField&,const GenericField&,const GenericField&);
    template<typename GenericField>
    void wenoloop3(lexer*,fdm*,const GenericField&,int,const GenericField&,const GenericField&,const GenericField&);
    template<typename GenericField>
    void wenoloop4(lexer*,fdm*,const GenericField&,int,const GenericField&,const GenericField&,const GenericField&);

    template<typename CF, typename MF>
    void aij_south(lexer*,fdm*,const CF&, MF&);
    template<typename CF, typename MF>
    void aij_north(lexer*,fdm*,const CF&, MF&);
    template<typename CF, typename MF>
    void aij_east(lexer*,fdm*,const CF&, MF&);
    template<typename CF, typename MF>
    void aij_west(lexer*,fdm*,const CF&, MF&);
    template<typename CF, typename MF>
    void aij_top(lexer*,fdm*,const CF&, MF&);
    template<typename CF, typename MF>
    void aij_bottom(lexer*,fdm*,const CF&, MF&);

    template<typename GenericField>
    void iqmin(lexer*, fdm*, const GenericField&);
    template<typename GenericField>
    void jqmin(lexer*, fdm*, const GenericField&);
    template<typename GenericField>
    void kqmin(lexer*, fdm*, const GenericField&);
    template<typename GenericField>
    void iqmax(lexer*, fdm*, const GenericField&);
    template<typename GenericField>
    void jqmax(lexer*, fdm*, const GenericField&);
    template<typename GenericField>
    void kqmax(lexer*, fdm*, const GenericField&);

    const double tttw,fourth,third,sevsix,elvsix,sixth,fivsix,tenth;
    const double sixten,treten,epsi,deltin;


    double is1,is2,is3;
    double alpha1,alpha2,alpha3;
    double w1,w2,w3;
    double umin, umax, uplus;
    int count;



    double ivel1,ivel2,jvel1,jvel2,kvel1,kvel2;
    double iadvec,jadvec,kadvec;

    flux *pflux;

    double *DX,*DY,*DZ;



};

#endif
