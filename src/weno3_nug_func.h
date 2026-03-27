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

#ifndef WENO3_NUG_FUNC_H_
#define WENO3_NUG_FUNC_H_

#include"increment.h"
#include"lexer.h"
#if USE_AMREX
#include <AMReX_GPU.H>
#endif

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

    // ------------------------------------------------------------------
    // Stateless face-flux kernels and divergence helpers.
    // All temporaries are local — no member state is read or written,
    // making these suitable for use in GPU device code.
    // fx_at / fy_at / fz_at compute the WENO flux at a single face.
    // fx_div / fy_div / fz_div combine both faces into the divergence
    // contribution vel_hi*f_hi - vel_lo*f_lo without mutating i/j/k.
    // ------------------------------------------------------------------
    #if USE_AMREX
        #define WENO3_NUG_GPU_ AMREX_GPU_HOST_DEVICE
    #else
        #define WENO3_NUG_GPU_
    #endif

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fx_at(GenericField& b, int i, int j, int k, double advec,
                                int ip, int uf, double eps, double ps)
    {
        if(advec>0.0)
        {
            const double q1=b(i-1,j,k), q2=b(i,j,k), q3=b(i+1,j,k);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfx[ip][uf][0]*(d1*d1), is2=isfx[ip][uf][1]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfx[ip][uf][0], c1=cfx[ip][uf][1], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfx[ip][uf][0][0]*q2+qfx[ip][uf][0][1]*q3)
                 + w2*(qfx[ip][uf][1][0]*q2-qfx[ip][uf][1][1]*q1);
        }
        else if(advec<0.0)
        {
            const double q1=b(i,j,k), q2=b(i+1,j,k), q3=b(i+2,j,k);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfx[ip][uf][2]*(d1*d1), is2=isfx[ip][uf][3]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfx[ip][uf][2], c1=cfx[ip][uf][3], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfx[ip][uf][2][0]*q2-qfx[ip][uf][2][1]*q3)
                 + w2*(qfx[ip][uf][3][0]*q1+qfx[ip][uf][3][1]*q2);
        }
        return 0.0;
    }

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fy_at(GenericField& b, int i, int j, int k, double advec,
                                int jp, int vf, double eps, double ps)
    {
        if(advec>0.0)
        {
            const double q1=b(i,j-1,k), q2=b(i,j,k), q3=b(i,j+1,k);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfy[jp][vf][0]*(d1*d1), is2=isfy[jp][vf][1]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfy[jp][vf][0], c1=cfy[jp][vf][1], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfy[jp][vf][0][0]*q2+qfy[jp][vf][0][1]*q3)
                 + w2*(qfy[jp][vf][1][0]*q2-qfy[jp][vf][1][1]*q1);
        }
        else if(advec<0.0)
        {
            const double q1=b(i,j,k), q2=b(i,j+1,k), q3=b(i,j+2,k);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfy[jp][vf][2]*(d1*d1), is2=isfy[jp][vf][3]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfy[jp][vf][2], c1=cfy[jp][vf][3], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfy[jp][vf][2][0]*q2-qfy[jp][vf][2][1]*q3)
                 + w2*(qfy[jp][vf][3][0]*q1+qfy[jp][vf][3][1]*q2);
        }
        return 0.0;
    }

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fz_at(GenericField& b, int i, int j, int k, double advec,
                                int kp, int wf, double eps, double ps)
    {
        if(advec>0.0)
        {
            const double q1=b(i,j,k-1), q2=b(i,j,k), q3=b(i,j,k+1);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfz[kp][wf][0]*(d1*d1), is2=isfz[kp][wf][1]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfz[kp][wf][0], c1=cfz[kp][wf][1], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfz[kp][wf][0][0]*q2+qfz[kp][wf][0][1]*q3)
                 + w2*(qfz[kp][wf][1][0]*q2-qfz[kp][wf][1][1]*q1);
        }
        else if(advec<0.0)
        {
            const double q1=b(i,j,k), q2=b(i,j,k+1), q3=b(i,j,k+2);
            const double d1=q3-q2, d2=q2-q1;
            const double is1=isfz[kp][wf][2]*(d1*d1), is2=isfz[kp][wf][3]*(d2*d2);
            const double s1sq=(is1+ps)*(is1+ps), s2sq=(is2+ps)*(is2+ps);
            const double c0=cfz[kp][wf][2], c1=cfz[kp][wf][3], a=c0/s1sq+c1/s2sq;
            const double w1=c0/(eps+s1sq*a), w2=c1/(eps+s2sq*a);
            return w1*(qfz[kp][wf][2][0]*q2-qfz[kp][wf][2][1]*q3)
                 + w2*(qfz[kp][wf][3][0]*q1+qfz[kp][wf][3][1]*q2);
        }
        return 0.0;
    }

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fx_div(GenericField& b, int i, int j, int k,
                                 double vel_lo, double vel_hi,
                                 int ip, int uf, double eps, double ps)
    {
        return vel_hi * fx_at(b, i,   j, k, vel_hi, ip,   uf, eps, ps)
             - vel_lo * fx_at(b, i-1, j, k, vel_lo, ip-1, uf, eps, ps);
    }

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fy_div(GenericField& b, int i, int j, int k,
                                 double vel_lo, double vel_hi,
                                 int jp, int vf, double eps, double ps)
    {
        return vel_hi * fy_at(b, i, j,   k, vel_hi, jp,   vf, eps, ps)
             - vel_lo * fy_at(b, i, j-1, k, vel_lo, jp-1, vf, eps, ps);
    }

    template<typename GenericField>
    WENO3_NUG_GPU_
    static inline double fz_div(GenericField& b, int i, int j, int k,
                                 double vel_lo, double vel_hi,
                                 int kp, int wf, double eps, double ps)
    {
        return vel_hi * fz_at(b, i, j, k,   vel_hi, kp,   wf, eps, ps)
             - vel_lo * fz_at(b, i, j, k-1, vel_lo, kp-1, wf, eps, ps);
    }

#undef WENO3_NUG_GPU_

    static double ****qfx,****qfy,****qfz;
    static double ***cfx,***cfy,***cfz;
    static double ***isfx,***isfy,***isfz;


    double q1,q2,q3;

    static constexpr double epsilon = 0.0, psi = 1.0e-6;
    double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
    double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;

    int uf,vf,wf;
private:

    lexer* p;

    static bool iniflag;
};

#endif
