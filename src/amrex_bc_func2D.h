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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef AMREX_BC_FUNC2D_H_
#define AMREX_BC_FUNC2D_H_

#include "definitions.h"

#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <initializer_list>
#include <type_traits>
#include <algorithm>

class amrex_bc_func2D
{
public:
    amrex_bc_func2D() = default;
    virtual ~amrex_bc_func2D() = default;

    enum class BoundaryConditionTypeLabel : int { NONE = 0, NEUMANN = 4, NOSLIP = 5, OUTFLOWBC = 6, SOMMERFELD = 7,
                                        POTENTIAL = 8, NEUMANN_X = 14, NEUMANN_HX = 41, NEUMANN_HY = 42};
    enum Gbc : int { INFLOW = 1, OUTFLOW = 2, SYMMETRY = 3, WAVEGEN = 6, NUMBEACH = 7, WALL = 21 };
private:
    enum Dir : int { X_NEG = 1, X_POS = 4, Y_NEG = 3, Y_POS = 2};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static bool is_cs_x(int cs)
    {
        return cs == Dir::X_NEG || cs == Dir::X_POS;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static bool is_cs_y(int cs)
    {
        return cs == Dir::Y_NEG || cs == Dir::Y_POS;
    }
public:
    struct Slice1BcDecision {
        struct Slice1Params {
            AMREX_GPU_HOST_DEVICE
            Slice1Params() noexcept = default;

            BoundaryConditionTypeLabel gclabel_u = BoundaryConditionTypeLabel::NEUMANN;
            bool awa_label = false;
            int B99 = 0;
        };

        AMREX_GPU_HOST_DEVICE
        Slice1BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Slice1BcDecision(const Slice1Params& params)
            : m_params(params)
        {}

        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            // general Neuman
            if(gcv==1 || gcv==40 || gcv==50)
                return BoundaryConditionTypeLabel::NEUMANN;
            //Wall
            // Parallel
            if((bc==NUMBEACH || bc==WALL) && is_cs_y(cs) && (gcv==1 || gcv==10 || gcv==20))
                return m_params.gclabel_u;
            // Orthogonal
            else if((bc==SYMMETRY || (bc==NUMBEACH && !m_params.awa_label) || bc==WALL) && is_cs_x(cs) && (gcv==1 || gcv==10 || gcv==20))
                return BoundaryConditionTypeLabel::NOSLIP;
            //Outflow
            else if(bc==OUTFLOW && is_cs_x(cs) && (gcv==1 || gcv==10 || gcv==20))
                return BoundaryConditionTypeLabel::OUTFLOWBC;
            //Symmetry
            else if(bc==SYMMETRY && is_cs_y(cs) && (gcv==1 || gcv==10 || gcv==20))
                return BoundaryConditionTypeLabel::NEUMANN;
            //Hx
            else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==52 || gcv==53))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==51 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN_HX;

            else if(bc==NUMBEACH && m_params.B99==3)
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;
            //Patch
            else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==1 || gcv==7 || gcv==10 || gcv==20))
                return BoundaryConditionTypeLabel::NEUMANN;
            //Patch Hx
            else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN_HX;

            else if((bc==222 || bc==212 || bc==122 || bc==112) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else
                return BoundaryConditionTypeLabel::NONE;
        };
    private:
        Slice1Params m_params{};
    };

    struct Slice2BcDecision {
        struct Slice2Params {
            AMREX_GPU_HOST_DEVICE
            Slice2Params() noexcept = default;

            BoundaryConditionTypeLabel gclabel_v = BoundaryConditionTypeLabel::NEUMANN;
            int B99 = 0;
        };

        AMREX_GPU_HOST_DEVICE
        Slice2BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Slice2BcDecision(const Slice2Params& params)
            : m_params(params)
        {}

        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            // general Neuman
            if(gcv==1 || gcv==40 || gcv==50)
                return BoundaryConditionTypeLabel::NEUMANN;

            //Wall
            // Parallel
            else if((bc==INFLOW || bc==WAVEGEN || bc==NUMBEACH || bc==WALL) && is_cs_x(cs) && (gcv==2 || gcv==11 || gcv==21))
                return m_params.gclabel_v;

            // Orthogonal
            else if((bc==SYMMETRY || bc==NUMBEACH || bc==WALL) && is_cs_y(cs) && (gcv==2 || gcv==11 || gcv==21))
                return BoundaryConditionTypeLabel::NOSLIP;

            //Patch
            else if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==2 || gcv==8 || gcv==11 || gcv==21))
                return BoundaryConditionTypeLabel::NEUMANN;

            //Outflow
            else if(bc==OUTFLOW && is_cs_y(cs) && (gcv==2 || gcv==11 || gcv==21))
                return BoundaryConditionTypeLabel::NEUMANN;

            // Symmetry
            else if(bc==SYMMETRY && is_cs_x(cs) && (gcv==2 || gcv==11 || gcv==21))
                return BoundaryConditionTypeLabel::NEUMANN;

            //Hy
            else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==51 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if(bc==NUMBEACH && (m_params.B99==3 || m_params.B99==4))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            //Patch Hy
            else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==55 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN_HY;

            else if((bc==222 || bc==212 || bc==122 || bc==112) && (gcv==55 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else
                return BoundaryConditionTypeLabel::NONE;

        };
    private:
        Slice2Params m_params{};
    };

    struct Slice4BcDecision {
        struct Slice4Params {
            AMREX_GPU_HOST_DEVICE
            Slice4Params() noexcept = default;

            int A515 = 0;
            int B98 = 0;
            int B99 = 0;
        };

        AMREX_GPU_HOST_DEVICE
        Slice4BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Slice4BcDecision(const Slice4Params& params)
            : m_params(params)
        {}

        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            // general Neuman
            if(gcv==1 || gcv==20 || gcv==30 || gcv==40 || gcv==50 || gcv==55)
                return BoundaryConditionTypeLabel::NEUMANN;

            // vertical w
            else if(bc!=INFLOW && (gcv==12 || gcv==24))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if(bc==INFLOW && gcv==12)
                return BoundaryConditionTypeLabel::NOSLIP;

            // pressure 40
            else if((bc==INFLOW || bc==OUTFLOW || bc==SYMMETRY || bc==WAVEGEN || bc==NUMBEACH || bc==WALL || bc==111 || bc==112 || bc==211 || bc==212) && gcv==44)
                return BoundaryConditionTypeLabel::NEUMANN;

            // Potential Ini
            else if((bc==SYMMETRY || bc==WALL) && gcv==49)
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==INFLOW || bc==OUTFLOW || bc==WAVEGEN || bc==NUMBEACH) && gcv==49)
                return BoundaryConditionTypeLabel::POTENTIAL;

            //Patch eta / Hx / Hy
            else if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==50 || gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            // eta
            else if((bc==INFLOW || bc==WAVEGEN) && (gcv==52 || gcv==54)  && (m_params.B98<3 || m_params.A515<=2))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if(bc==OUTFLOW && (gcv==51 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if((bc==SYMMETRY || bc==WALL) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if(bc==NUMBEACH && (gcv==51 || gcv==54 || m_params.B99==0 || m_params.B99==1 || ((gcv==52 || gcv==53) && m_params.B99==3)))
                return BoundaryConditionTypeLabel::NEUMANN;

            else if(bc==NUMBEACH && (gcv==51 || gcv==52 || gcv==53 || gcv==54) && m_params.B99==4)
                return BoundaryConditionTypeLabel::SOMMERFELD;

            // Fifsf 60 - 3D
            else if(((cs==X_NEG && m_params.B98<=2) || cs==Y_POS || cs==Y_NEG || (cs==X_POS && m_params.B99<=2)) && gcv==60)
                return BoundaryConditionTypeLabel::NEUMANN;

            // eta 150
            else if((bc==INFLOW || bc==WAVEGEN) && (gcv==152 || gcv==154))
                return BoundaryConditionTypeLabel::NEUMANN_X;

            else if((bc==OUTFLOW || bc==NUMBEACH) && (gcv==151 || gcv==154))
                return BoundaryConditionTypeLabel::NEUMANN_X;

            else if(bc==NUMBEACH && (gcv==151 || gcv==152 || gcv==153 || gcv==154) &&m_params.B99==3)
                return BoundaryConditionTypeLabel::NEUMANN_X;

            else if(bc==NUMBEACH && (gcv==151 || gcv==152 || gcv==153 || gcv==154) &&m_params.B99==4)
                return BoundaryConditionTypeLabel::SOMMERFELD;

            else if((bc==SYMMETRY || bc==WALL) && (gcv==151 || gcv==152 || gcv==153 || gcv==154))
                return BoundaryConditionTypeLabel::NEUMANN_X;

            else if(gcv==155)
                return BoundaryConditionTypeLabel::NEUMANN_X;

            // Fifsf 160 - 2D
            else if(((cs==X_NEG && m_params.B98<=2) || (cs==X_POS && m_params.B99<=2)) && gcv==160)
                return BoundaryConditionTypeLabel::NEUMANN_X;

            else
                return BoundaryConditionTypeLabel::NONE;
        };
    private:
        Slice4Params m_params{};
    };

    struct ConstMyExtBCFillSliceParams {
        AMREX_GPU_HOST_DEVICE
        ConstMyExtBCFillSliceParams() noexcept = default;
        ConstMyExtBCFillSliceParams(const amrex::Array<int,4>& bc_values_in,
                                    bool y_dimension_exists_in, DataLocation data_location_in, double gravity_in) noexcept
            : bc_values(bc_values_in),
              y_dimension_exists(y_dimension_exists_in), data_location(data_location_in), gravity(gravity_in) {}

        // Per-side Gbc values used by the BC decision.
        // Indexed as: 0=X_NEG  1=X_POS  2=Y_NEG  3=Y_POS  (same order as face_labels)
        const amrex::Array<int,4> bc_values = {};
        bool y_dimension_exists = true;
        const DataLocation data_location = DataLocation::CELL_CENTERED;
        double gravity = 9.81;
    };

    struct MyExtBCFillSliceParams {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillSliceParams() noexcept = default;
        MyExtBCFillSliceParams(double Ui_in, double Uo_in, double dt_in, double wd_in, const amrex::Array<int,4>& face_labels_in) noexcept
            : Ui(Ui_in), Uo(Uo_in), dt(dt_in), wd(wd_in), face_labels(face_labels_in) {}

        double Ui = 0.0;
        double Uo = 0.0;
        double dt = 0.0;
        double wd = 0.0;
        // Path B: pre-evaluated face labels for single-component fields.
        // Indexed as: 0=X_NEG  1=X_POS  2=Y_NEG  3=Y_POS
        // Populated in FillDomainBoundaryImpl; avoids BCRec pointer dereference
        // in the common (numcomp==1) face path.
        amrex::Array<int,4> face_labels = {};
    };

    template <typename BCDecision>
    struct MyExtBCFillSlice
    {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillSlice() = delete;

        AMREX_GPU_HOST_DEVICE
        MyExtBCFillSlice(const ConstMyExtBCFillSliceParams const_params, const MyExtBCFillSliceParams params)
            : m_const_params(const_params), m_params(params) {}

        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                        const int dcomp, const int numcomp,
                        amrex::GeometryData const& geom, const amrex::Real time,
                        const amrex::BCRec* bcr, const int bcomp,
                        const int orig_comp) const
        {
            // slice_amrex is constructed without ghost cells in z-dir,
            // so this should not needed
            if(iv[2]!=0)
            {
                for(int n=0; n<numcomp; ++n)
                    dest(iv, dcomp+n) = amrex::Real(0);
                return;
            }

            if(!m_const_params.y_dimension_exists && iv[1]!=0)
            {
                for(int n=0; n<numcomp; ++n)
                    dest(iv, dcomp+n) = amrex::Real(0);
                return;
            }

            const amrex::Box dom = geom.Domain();

            // Path D: single pass computes both the interior clamp and the
            // out-of-bounds flags, halving the comparisons vs. two separate passes.
            amrex::IntVect interior = iv;
            int out_x = 0, out_y = 0;

            if (iv[0] < dom.smallEnd(0)) {
                out_x = -1; interior[0] = dom.smallEnd(0);
            } else if (iv[0] > dom.bigEnd(0)) {
                out_x = 1; interior[0] = dom.bigEnd(0);
                if (m_const_params.data_location == DataLocation::FACE_X) --interior[0];
            }

            if (iv[1] < dom.smallEnd(1)) {
                out_y = -1; interior[1] = dom.smallEnd(1);
            } else if (iv[1] > dom.bigEnd(1)) {
                out_y = 1; interior[1] = dom.bigEnd(1);
                if (m_const_params.data_location == DataLocation::FACE_Y) --interior[1];
            }

            if (out_x == 0 || out_y == 0)
            {
                // Pure face cell.  Encode face as an index matching face_labels:
                //   0=X_NEG  1=X_POS  2=Y_NEG  3=Y_POS
                int face_idx, cs;
                if      (out_x == -1) { face_idx = 0; cs = Dir::X_NEG; }
                else if (out_x ==  1) { face_idx = 1; cs = Dir::X_POS; }
                else if (out_y == -1) { face_idx = 2; cs = Dir::Y_NEG; }
                else                  { face_idx = 3; cs = Dir::Y_POS; }

                // Path B: single-component fields (the common case) read their
                // label from the inline face_labels array stored in the functor,
                // avoiding a BCRec pointer dereference.  Multi-component batch
                // calls derive dim/side from face_idx and use the BCRec array.
                if (numcomp == 1)
                {
                    const BoundaryConditionTypeLabel label = static_cast<BoundaryConditionTypeLabel>(
                        m_params.face_labels[face_idx]);
                    apply_label(label, cs, dest, iv, interior, dcomp, geom);
                }
                else
                {
                    const int  bcr_dim = (face_idx < 2) ? 0 : 1;
                    const bool bcr_lo  = (face_idx % 2 == 0);
                    for(int n = 0; n < numcomp; ++n)
                    {
                        const BoundaryConditionTypeLabel label = static_cast<BoundaryConditionTypeLabel>(
                            bcr_lo ? bcr[bcomp+n].lo(bcr_dim) : bcr[bcomp+n].hi(bcr_dim));
                        apply_label(label, cs, dest, iv, interior, dcomp+n, geom);
                    }
                }
                return;
            }

            // X+Y Corners are the only remaining out-of-bounds case.
            // They are either Neumann in layer 1 or NOSLIP in layer 2.
            for (int n = 0; n < numcomp; ++n)
            {
                BoundaryConditionTypeLabel label = is_corner_layer1(iv, dom)
                        ? BoundaryConditionTypeLabel::NEUMANN
                        : BoundaryConditionTypeLabel::NOSLIP;

                apply_label(label, 0, dest, iv, interior, dcomp+n, geom);
            }
        };
    private:
        // Apply a single BC label to one component of a ghost cell.
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void apply_label(BoundaryConditionTypeLabel label, int cs,
                         amrex::Array4<amrex::Real> const& dest,
                         const amrex::IntVect& iv, const amrex::IntVect& interior,
                         int comp, amrex::GeometryData const& geom) const
        {
            switch (label)
            {
                case BoundaryConditionTypeLabel::NONE:
                    break;
                case BoundaryConditionTypeLabel::NEUMANN:
                default:
                    dest(iv, comp) = dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::NEUMANN_X:
                    if(is_cs_x(cs))
                        dest(iv, comp) = dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::NEUMANN_HX:
                    if(is_cs_x(cs))
                    {
                        if(abs(iv[0]-interior[0])>1)
                            dest(iv, comp) = dest(interior, comp);
                    }
                    else
                        dest(iv, comp) = dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::NEUMANN_HY:
                    if(is_cs_y(cs))
                    {
                        if(abs(iv[1]-interior[1])>1)
                            dest(iv, comp) = dest(interior, comp);
                    }
                    else
                        dest(iv, comp) = dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::NOSLIP:
                    dest(iv, comp) = amrex::Real(0);
                    break;
                case BoundaryConditionTypeLabel::OUTFLOWBC: // only possible for x-dir
                    if(cs==Dir::X_POS)
                        dest(iv, comp) = std::max(amrex::Real(0), dest(interior, comp));
                    else
                        dest(iv, comp) = dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::POTENTIAL:
                    if(cs==Dir::X_NEG)
                        dest(iv, comp) = dest(interior, comp) - m_params.Ui * geom.CellSize(0);
                    else if(cs==Dir::X_POS)
                        dest(iv, comp) = dest(interior, comp) + m_params.Uo * geom.CellSize(0);
                    break;

                case BoundaryConditionTypeLabel::SOMMERFELD:
                    if(cs==Dir::X_NEG)
                    {
                        const double interior_value = dest(interior, comp);
                        dest(iv, comp) = interior_value - m_params.dt * sqrt(m_const_params.gravity * (m_params.wd + interior_value)) * (dest(interior+amrex::IntVect(1,0,0), comp) - interior_value)/geom.CellSize(0);
                    }
                    else if(cs==Dir::X_POS)
                    {
                        const double interior_value = dest(interior, comp);
                        dest(iv, comp) = interior_value - m_params.dt * sqrt(m_const_params.gravity * (m_params.wd + interior_value)) * (dest(interior+amrex::IntVect(-1,0,0), comp) - interior_value)/geom.CellSize(0);
                    }
                    else if(cs==Dir::Y_NEG)
                    {
                        const double interior_value = dest(interior, comp);
                        dest(iv, comp) = interior_value - m_params.dt * sqrt(m_const_params.gravity * (m_params.wd + interior_value)) * (dest(interior+amrex::IntVect(0,1,0), comp) - interior_value)/geom.CellSize(0);
                    }
                    else if(cs==Dir::Y_POS)
                    {
                        const double interior_value = dest(interior, comp);
                        dest(iv, comp) = interior_value - m_params.dt * sqrt(m_const_params.gravity * (m_params.wd + interior_value)) * (dest(interior+amrex::IntVect(0,-1,0), comp) - interior_value)/geom.CellSize(0);
                    }
                    break;
            }
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_corner_layer1(const amrex::IntVect& iv, const amrex::Box& box) const
        {
            int max_dist = 0;
            for(int dir=0; dir<2; ++dir)
            {
                int d = 0;
                if (iv[dir] < box.smallEnd(dir))
                {
                    d = box.smallEnd(dir) - iv[dir];
                }
                else if (iv[dir] > box.bigEnd(dir))
                {
                    d = iv[dir] - box.bigEnd(dir);
                }
                if (d > max_dist) max_dist = d;
            }

            return (max_dist == 1);
        }

        const ConstMyExtBCFillSliceParams m_const_params{};
        const MyExtBCFillSliceParams m_params{};
    };
};
#endif
#endif
