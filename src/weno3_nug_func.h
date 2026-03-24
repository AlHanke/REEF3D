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

#ifndef WENO3_NUG_FUNC_H_
#define WENO3_NUG_FUNC_H_

#include"increment.h"
#include"lexer.h"

class weno3_nug_func : public increment
{
public:
    weno3_nug_func(lexer*);
    virtual ~weno3_nug_func() = default;

    void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);

    void ini(lexer*);

    inline void is_min_x()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1x=isfx[IP][uf][0]*(d1*d1);
        is2x=isfx[IP][uf][1]*(d2*d2);
    }
    inline void is_max_x()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1x=isfx[IP][uf][2]*(d1*d1);
        is2x=isfx[IP][uf][3]*(d2*d2);
    }
    inline void is_min_y()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1y=isfy[JP][vf][0]*(d1*d1);
        is2y=isfy[JP][vf][1]*(d2*d2);
    }
    inline void is_max_y()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1y=isfy[JP][vf][2]*(d1*d1);
        is2y=isfy[JP][vf][3]*(d2*d2);
    }
    inline void is_min_z()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1z=isfz[KP][wf][0]*(d1*d1);
        is2z=isfz[KP][wf][1]*(d2*d2);
    }
    inline void is_max_z()
    {
        const double d1=q3-q2, d2=q2-q1;
        is1z=isfz[KP][wf][2]*(d1*d1);
        is2z=isfz[KP][wf][3]*(d2*d2);
    }

    inline void weight_min_x()
    {
        const double s1sq=(is1x+psi)*(is1x+psi), s2sq=(is2x+psi)*(is2x+psi);
        const double c0=cfx[IP][uf][0], c1=cfx[IP][uf][1], a=c0/s1sq+c1/s2sq;
        w1x=c0/(epsilon+s1sq*a); w2x=c1/(epsilon+s2sq*a);
    }
    inline void weight_max_x()
    {
        const double s1sq=(is1x+psi)*(is1x+psi), s2sq=(is2x+psi)*(is2x+psi);
        const double c0=cfx[IP][uf][2], c1=cfx[IP][uf][3], a=c0/s1sq+c1/s2sq;
        w1x=c0/(epsilon+s1sq*a); w2x=c1/(epsilon+s2sq*a);
    }
    inline void weight_min_y()
    {
        const double s1sq=(is1y+psi)*(is1y+psi), s2sq=(is2y+psi)*(is2y+psi);
        const double c0=cfy[JP][vf][0], c1=cfy[JP][vf][1], a=c0/s1sq+c1/s2sq;
        w1y=c0/(epsilon+s1sq*a); w2y=c1/(epsilon+s2sq*a);
    }
    inline void weight_max_y()
    {
        const double s1sq=(is1y+psi)*(is1y+psi), s2sq=(is2y+psi)*(is2y+psi);
        const double c0=cfy[JP][vf][2], c1=cfy[JP][vf][3], a=c0/s1sq+c1/s2sq;
        w1y=c0/(epsilon+s1sq*a); w2y=c1/(epsilon+s2sq*a);
    }
    inline void weight_min_z()
    {
        const double s1sq=(is1z+psi)*(is1z+psi), s2sq=(is2z+psi)*(is2z+psi);
        const double c0=cfz[KP][wf][0], c1=cfz[KP][wf][1], a=c0/s1sq+c1/s2sq;
        w1z=c0/(epsilon+s1sq*a); w2z=c1/(epsilon+s2sq*a);
    }
    inline void weight_max_z()
    {
        const double s1sq=(is1z+psi)*(is1z+psi), s2sq=(is2z+psi)*(is2z+psi);
        const double c0=cfz[KP][wf][2], c1=cfz[KP][wf][3], a=c0/s1sq+c1/s2sq;
        w1z=c0/(epsilon+s1sq*a); w2z=c1/(epsilon+s2sq*a);
    }

    static double ****qfx,****qfy,****qfz;
    static double ***cfx,***cfy,***cfz;
    static double ***isfx,***isfy,***isfz;

    static bool iniflag;

    double q1,q2,q3;

    const double epsilon,psi;
    double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
    double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;

    int uf,vf,wf;
private:

    lexer* p;
};

#endif
