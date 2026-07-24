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

#ifndef WENO_NUG_FUNC_H_
#define WENO_NUG_FUNC_H_

#if defined(_MSC_VER)
    #define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FORCE_INLINE __attribute__((always_inline)) inline
#else
    #define FORCE_INLINE inline
#endif

#include "increment.h"
#include "lexer.h"
#include "field.h"
#include "slice.h"
#include <array>
#include <vector>

using namespace std;

class weno_nug_func : public increment
{
public:
    weno_nug_func(lexer*);
    virtual ~weno_nug_func();

    void precalc_qf(lexer*);
    void precalc_cf(lexer*);
    void precalc_isf(lexer*);

    void ini(lexer*);
    void rebuild_levels(lexer* p, int new_nlevs);

    // IS ----
    // x
    FORCE_INLINE void is_min_x() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1x = isfx[IP][uf][0][0]*dq54*dq54 + isfx[IP][uf][0][1]*(dq54)*(dq34) + isfx[IP][uf][0][2]*dq34*dq34;
        is2x = isfx[IP][uf][1][0]*dq23*dq23 + isfx[IP][uf][1][1]*(dq43)*(dq23) + isfx[IP][uf][1][2]*dq43*dq43;
        is3x = isfx[IP][uf][2][0]*dq12*dq12 + isfx[IP][uf][2][1]*(dq32)*(dq12) + isfx[IP][uf][2][2]*dq32*dq32;
    }
    FORCE_INLINE void is_max_x() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1x = isfx[IP][uf][3][0]*dq54*dq54 + isfx[IP][uf][3][1]*(dq54)*(dq34) + isfx[IP][uf][3][2]*dq34*dq34;
        is2x = isfx[IP][uf][4][0]*dq23*dq23 + isfx[IP][uf][4][1]*(dq43)*(dq23) + isfx[IP][uf][4][2]*dq43*dq43;
        is3x = isfx[IP][uf][5][0]*dq12*dq12 + isfx[IP][uf][5][1]*(dq32)*(dq12) + isfx[IP][uf][5][2]*dq32*dq32;
    }

    // y
    FORCE_INLINE void is_min_y() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1y = isfy[JP][vf][0][0]*dq54*dq54 + isfy[JP][vf][0][1]*(dq54)*(dq34) + isfy[JP][vf][0][2]*dq34*dq34;
        is2y = isfy[JP][vf][1][0]*dq23*dq23 + isfy[JP][vf][1][1]*(dq43)*(dq23) + isfy[JP][vf][1][2]*dq43*dq43;
        is3y = isfy[JP][vf][2][0]*dq12*dq12 + isfy[JP][vf][2][1]*(dq32)*(dq12) + isfy[JP][vf][2][2]*dq32*dq32;
    }
    FORCE_INLINE void is_max_y() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1y = isfy[JP][vf][3][0]*dq54*dq54 + isfy[JP][vf][3][1]*(dq54)*(dq34) + isfy[JP][vf][3][2]*dq34*dq34;
        is2y = isfy[JP][vf][4][0]*dq23*dq23 + isfy[JP][vf][4][1]*(dq43)*(dq23) + isfy[JP][vf][4][2]*dq43*dq43;
        is3y = isfy[JP][vf][5][0]*dq12*dq12 + isfy[JP][vf][5][1]*(dq32)*(dq12) + isfy[JP][vf][5][2]*dq32*dq32;
    }

    // z
    FORCE_INLINE void is_min_z() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1z = isfz[KP][wf][0][0]*dq54*dq54 + isfz[KP][wf][0][1]*(dq54)*(dq34) + isfz[KP][wf][0][2]*dq34*dq34;
        is2z = isfz[KP][wf][1][0]*dq23*dq23 + isfz[KP][wf][1][1]*(dq43)*(dq23) + isfz[KP][wf][1][2]*dq43*dq43;
        is3z = isfz[KP][wf][2][0]*dq12*dq12 + isfz[KP][wf][2][1]*(dq32)*(dq12) + isfz[KP][wf][2][2]*dq32*dq32;
    }
    FORCE_INLINE void is_max_z() noexcept
    {
        const double dq12 = q1 - q2;
        const double dq23 = q2 - q3;
        const double dq32 = q3 - q2;
        const double dq34 = q3 - q4;
        const double dq43 = q4 - q3;
        const double dq54 = q5 - q4;

        is1z = isfz[KP][wf][3][0]*dq54*dq54 + isfz[KP][wf][3][1]*(dq54)*(dq34) + isfz[KP][wf][3][2]*dq34*dq34;
        is2z = isfz[KP][wf][4][0]*dq23*dq23 + isfz[KP][wf][4][1]*(dq43)*(dq23) + isfz[KP][wf][4][2]*dq43*dq43;
        is3z = isfz[KP][wf][5][0]*dq12*dq12 + isfz[KP][wf][5][1]*(dq32)*(dq12) + isfz[KP][wf][5][2]*dq32*dq32;
    }

    // Weights ----
    //
    // The WENO-JS weights are w_i = (c_i/d_i) / sum_j(c_j/d_j), with d_i = (is_i+psi)^2.
    // Written that way each weight costs a division, plus one more for the normalisation.
    // Multiplying numerator and denominator through by d1*d2*d3 gives the algebraically
    // identical form
    //     n_i = c_i * d_j * d_k        (j,k the other two)
    //     w_i = n_i / (n1 + n2 + n3)
    // which needs a single division. Trades 3 divisions for 6 multiplications per call;
    // at 6 calls per cell that is 24 divisions/cell down to 6.
    //
    // Range: d_i ~ (dq)^4, so n_i ~ (dq)^8. Overflow needs |dq| ~ 1e38, underflow to zero
    // needs |dq| ~ 1e-40; neither is reachable for physical field differences. The ratio is
    // exact in real arithmetic but rounds differently from the original, so results are
    // close to but not bit-identical with the previous form.
    // x
    FORCE_INLINE void weight_min_x() noexcept
    {
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;

        const double d1x = is1x_psi*is1x_psi;
        const double d2x = is2x_psi*is2x_psi;
        const double d3x = is3x_psi*is3x_psi;

        const double n1x = cfx[IP][uf][0]*d2x*d3x;
        const double n2x = cfx[IP][uf][1]*d1x*d3x;
        const double n3x = cfx[IP][uf][2]*d1x*d2x;

        const double inv_sumx = 1.0/(n1x + n2x + n3x);
        w1x = n1x*inv_sumx;
        w2x = n2x*inv_sumx;
        w3x = n3x*inv_sumx;
    }
    FORCE_INLINE void weight_max_x() noexcept
    {
        const double is1x_psi = is1x + psi;
        const double is2x_psi = is2x + psi;
        const double is3x_psi = is3x + psi;

        const double d1x = is1x_psi*is1x_psi;
        const double d2x = is2x_psi*is2x_psi;
        const double d3x = is3x_psi*is3x_psi;

        const double n1x = cfx[IP][uf][3]*d2x*d3x;
        const double n2x = cfx[IP][uf][4]*d1x*d3x;
        const double n3x = cfx[IP][uf][5]*d1x*d2x;

        const double inv_sumx = 1.0/(n1x + n2x + n3x);
        w1x = n1x*inv_sumx;
        w2x = n2x*inv_sumx;
        w3x = n3x*inv_sumx;
    }

    // y
    FORCE_INLINE void weight_min_y() noexcept
    {
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;

        const double d1y = is1y_psi*is1y_psi;
        const double d2y = is2y_psi*is2y_psi;
        const double d3y = is3y_psi*is3y_psi;

        const double n1y = cfy[JP][vf][0]*d2y*d3y;
        const double n2y = cfy[JP][vf][1]*d1y*d3y;
        const double n3y = cfy[JP][vf][2]*d1y*d2y;

        const double inv_sumy = 1.0/(n1y + n2y + n3y);
        w1y = n1y*inv_sumy;
        w2y = n2y*inv_sumy;
        w3y = n3y*inv_sumy;
    }
    FORCE_INLINE void weight_max_y() noexcept
    {
        const double is1y_psi = is1y + psi;
        const double is2y_psi = is2y + psi;
        const double is3y_psi = is3y + psi;

        const double d1y = is1y_psi*is1y_psi;
        const double d2y = is2y_psi*is2y_psi;
        const double d3y = is3y_psi*is3y_psi;

        const double n1y = cfy[JP][vf][3]*d2y*d3y;
        const double n2y = cfy[JP][vf][4]*d1y*d3y;
        const double n3y = cfy[JP][vf][5]*d1y*d2y;

        const double inv_sumy = 1.0/(n1y + n2y + n3y);
        w1y = n1y*inv_sumy;
        w2y = n2y*inv_sumy;
        w3y = n3y*inv_sumy;
    }

    // z
    FORCE_INLINE void weight_min_z() noexcept
    {
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;

        const double d1z = is1z_psi*is1z_psi;
        const double d2z = is2z_psi*is2z_psi;
        const double d3z = is3z_psi*is3z_psi;

        const double n1z = cfz[KP][wf][0]*d2z*d3z;
        const double n2z = cfz[KP][wf][1]*d1z*d3z;
        const double n3z = cfz[KP][wf][2]*d1z*d2z;

        const double inv_sumz = 1.0/(n1z + n2z + n3z);
        w1z = n1z*inv_sumz;
        w2z = n2z*inv_sumz;
        w3z = n3z*inv_sumz;
    }
    FORCE_INLINE void weight_max_z() noexcept
    {
        const double is1z_psi = is1z + psi;
        const double is2z_psi = is2z + psi;
        const double is3z_psi = is3z + psi;

        const double d1z = is1z_psi*is1z_psi;
        const double d2z = is2z_psi*is2z_psi;
        const double d3z = is3z_psi*is3z_psi;

        const double n1z = cfz[KP][wf][3]*d2z*d3z;
        const double n2z = cfz[KP][wf][4]*d1z*d3z;
        const double n3z = cfz[KP][wf][5]*d1z*d2z;

        const double inv_sumz = 1.0/(n1z + n2z + n3z);
        w1z = n1z*inv_sumz;
        w2z = n2z*inv_sumz;
        w3z = n3z*inv_sumz;
    }

    static inline std::vector<std::array<std::array<std::array<double, 2>, 6>, 2>> qfx, qfy, qfz;
    static inline std::vector<std::array<std::array<double, 6>, 2>> cfx, cfy, cfz;
    static inline std::vector<std::array<std::array<std::array<double, 3>, 6>, 2>> isfx, isfy, isfz;

    double q1,q2,q3,q4,q5;

    double is1x,is2x,is3x;
    double is1y,is2y,is3y;
    double is1z,is2z,is3z;
    double w1x,w2x,w3x;
    double w1y,w2y,w3y;
    double w1z,w2z,w3z;

    int uf,vf,wf;
protected:
    template<typename GenericField> FORCE_INLINE void iqmin(const GenericField& f) noexcept
    {
        const double v0=f(i-3,j,k), v1=f(i-2,j,k), v2=f(i-1,j,k),
                     v3=f(i,j,k),   v4=f(i+1,j,k), v5=f(i+2,j,k);
        q1 = (v1-v0)/p->DXP[IM3];
        q2 = (v2-v1)/p->DXP[IM2];
        q3 = (v3-v2)/p->DXP[IM1];
        q4 = (v4-v3)/p->DXP[IP];
        q5 = (v5-v4)/p->DXP[IP1];
    }
    template<typename GenericField> FORCE_INLINE void iqmax(const GenericField& f) noexcept
    {
        const double v0=f(i-2,j,k), v1=f(i-1,j,k), v2=f(i,j,k),
                     v3=f(i+1,j,k), v4=f(i+2,j,k), v5=f(i+3,j,k);
        q1 = (v1-v0)/p->DXP[IM2];
        q2 = (v2-v1)/p->DXP[IM1];
        q3 = (v3-v2)/p->DXP[IP];
        q4 = (v4-v3)/p->DXP[IP1];
        q5 = (v5-v4)/p->DXP[IP2];
    }

    template<typename GenericField> FORCE_INLINE void jqmin(const GenericField& f) noexcept
    {
        const double v0=f(i,j-3,k), v1=f(i,j-2,k), v2=f(i,j-1,k),
                     v3=f(i,j,k),   v4=f(i,j+1,k), v5=f(i,j+2,k);
        q1 = (v1-v0)/p->DYP[JM3];
        q2 = (v2-v1)/p->DYP[JM2];
        q3 = (v3-v2)/p->DYP[JM1];
        q4 = (v4-v3)/p->DYP[JP];
        q5 = (v5-v4)/p->DYP[JP1];
    }
    template<typename GenericField> FORCE_INLINE void jqmax(const GenericField& f) noexcept
    {
        const double v0=f(i,j-2,k), v1=f(i,j-1,k), v2=f(i,j,k),
                     v3=f(i,j+1,k), v4=f(i,j+2,k), v5=f(i,j+3,k);
        q1 = (v1-v0)/p->DYP[JM2];
        q2 = (v2-v1)/p->DYP[JM1];
        q3 = (v3-v2)/p->DYP[JP];
        q4 = (v4-v3)/p->DYP[JP1];
        q5 = (v5-v4)/p->DYP[JP2];
    }

    template<typename GenericField> FORCE_INLINE void kqmin(const GenericField& f) noexcept
    {
        const double v0=f(i,j,k-3), v1=f(i,j,k-2), v2=f(i,j,k-1),
                     v3=f(i,j,k),   v4=f(i,j,k+1), v5=f(i,j,k+2);
        q1 = (v1-v0)/p->DZP[KM3];
        q2 = (v2-v1)/p->DZP[KM2];
        q3 = (v3-v2)/p->DZP[KM1];
        q4 = (v4-v3)/p->DZP[KP];
        q5 = (v5-v4)/p->DZP[KP1];
    }
    template<typename GenericField> FORCE_INLINE void kqmax(const GenericField& f) noexcept
    {
        const double v0=f(i,j,k-2), v1=f(i,j,k-1), v2=f(i,j,k),
                     v3=f(i,j,k+1), v4=f(i,j,k+2), v5=f(i,j,k+3);
        q1 = (v1-v0)/p->DZP[KM2];
        q2 = (v2-v1)/p->DZP[KM1];
        q3 = (v3-v2)/p->DZP[KP];
        q4 = (v4-v3)/p->DZP[KP1];
        q5 = (v5-v4)/p->DZP[KP2];
    }

    FORCE_INLINE void isqmin(slice& f)
    {
        q1 = (f(i-2,j)-f(i-3,j))/p->DXP[IM3];
        q2 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
        q3 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
        q4 = (f(i+1,j)-f(i,j))/p->DXP[IP];
        q5 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
    }
    FORCE_INLINE void isqmax(slice& f)
    {
        q1 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];
        q2 = (f(i,j)-f(i-1,j))/p->DXP[IM1];
        q3 = (f(i+1,j)-f(i,j))/p->DXP[IP];
        q4 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
        q5 = (f(i+3,j)-f(i+2,j))/p->DXP[IP2];
    }

    FORCE_INLINE void jsqmin(slice& f)
    {
        q1 = (f(i,j-2)-f(i,j-3))/p->DYP[JM3];
        q2 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
        q3 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
        q4 = (f(i,j+1)-f(i,j))/p->DYP[JP];
        q5 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
    }
    FORCE_INLINE void jsqmax(slice& f)
    {
        q1 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];
        q2 = (f(i,j)-f(i,j-1))/p->DYP[JM1];
        q3 = (f(i,j+1)-f(i,j))/p->DYP[JP];
        q4 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
        q5 = (f(i,j+3)-f(i,j+2))/p->DYP[JP2];
    }

    static constexpr double psi = 1.0e-6;
private:
    static inline bool iniflag = false;

    int i_size, j_size, k_size;
    lexer *p;

    #if USE_AMREX
    friend class grid_amrex;
    #endif
};

#undef FORCE_INLINE

#endif
