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

#ifndef WENO_FLUX_NUG_H_
#define WENO_FLUX_NUG_H_

#include"convection.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"
#include"weno_nug_func.h"

// Defined after the includes: weno_nug_func.h #undef's FORCE_INLINE at its end,
// which would otherwise remove it again before the declarations below use it.
#if defined(_MSC_VER)
    #define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
    #define FORCE_INLINE __attribute__((always_inline)) inline
#else
    #define FORCE_INLINE inline
#endif

class weno_flux_nug final : public convection, public weno_nug_func
{
public:
    weno_flux_nug(lexer*);
    virtual ~weno_flux_nug() = default;

    void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    template<typename FluxT, typename GenericField>
    inline double aij(FluxT&, lexer*, fdm*, const GenericField&, int, const GenericField&, const GenericField&, const GenericField&, double*, double*, double*);

    template<typename GenericField>
    FORCE_INLINE double fx(lexer*, fdm*, const GenericField&, const GenericField&, double, int di=0);
    template<typename GenericField>
    FORCE_INLINE double fy(lexer*, fdm*, const GenericField&, const GenericField&, double, int dj=0);
    template<typename GenericField>
    FORCE_INLINE double fz(lexer*, fdm*, const GenericField&, const GenericField&, double, int dk=0);

    std::variant<flux_face_CDS2, flux_face_FOU,
                 flux_face_CDS2_vrans, flux_face_FOU_vrans,
                 flux_face_CDS2_2D, flux_face_FOU_2D,
                 flux_face_CDS2_vrans_2D, flux_face_FOU_vrans_2D> pflux;
};

#undef FORCE_INLINE

#endif
