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
#ifndef AMREX_BC_FUNC_H_
#define AMREX_BC_FUNC_H_

#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <initializer_list>
#include <type_traits>
#include <algorithm>

class amrex_bc_func
{
public:
    amrex_bc_func() = default;
    virtual ~amrex_bc_func() = default;

    enum class BoundaryConditionTypeLabel : int { NONE = 0, DIRICHLET_ORTH = 1, NEUMANN = 4, NOSLIP = 5, OUTFLOWBC = 6, SOMMERFELD = 7,
                                        POTENTIAL = 8, DIRICHLET_ORTH_REFLECT = 11, DIRICHLET_PARA_REFLECT = 12,
                                        NEUMANN_X = 14, NEUMANN_HX = 41, NEUMANN_HY = 42, HEATBC = 61 };
    enum class DataLocation : unsigned int { CELL_CENTERED = 0, FACE_X = 1, FACE_Y = 2, FACE_Z = 3 };
    enum Gbc : int { INFLOW = 1, OUTFLOW = 2, SYMMETRY = 3, WAVEGEN = 6, NUMBEACH = 7, WALL = 21 };
private:
    enum Dir : int { X_NEG = 1, X_POS = 4, Y_NEG = 3, Y_POS = 2, Z_NEG = 5, Z_POS = 6 };
public:
    struct Field1BcDecision {
        struct Field1Params {
            AMREX_GPU_HOST_DEVICE
            Field1Params() noexcept = default;

            BoundaryConditionTypeLabel gclabel_u = BoundaryConditionTypeLabel::NEUMANN;
            BoundaryConditionTypeLabel orth_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel inflow_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel outflow_label = BoundaryConditionTypeLabel::NEUMANN;
            bool awa_label = false;
            bool gclabel_outflow = true;
            bool i10_enabled = false;
        };

        AMREX_GPU_HOST_DEVICE
        Field1BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Field1BcDecision(const Field1Params& params)
            : m_params(params)
        {}

        AMREX_GPU_DEVICE
        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            if (gcv == 50)
                return BoundaryConditionTypeLabel::NEUMANN;

            const bool is_parallel_wall = ((bc == Gbc::NUMBEACH && !m_params.awa_label) || bc == Gbc::WALL);
            if (is_parallel_wall && is_cs_yz(cs) && matches_gcv(gcv, {1, 10, 114}))
                return m_params.gclabel_u;

            if (is_parallel_wall && is_cs_yz(cs) && gcv == 110)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == Gbc::WALL && gcv == 14)
                return BoundaryConditionTypeLabel::NEUMANN;

            if ((bc == Gbc::SYMMETRY || bc == Gbc::WALL) && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return m_params.orth_label;

            if (bc == Gbc::WAVEGEN && is_cs_x(cs) && matches_gcv(gcv, {1, 7, 10}))
                return m_params.inflow_label;

            if (bc == Gbc::WALL && is_cs_x(cs) && gcv == 7)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {1, 7, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (((bc == Gbc::OUTFLOW && m_params.gclabel_outflow) || bc == Gbc::SYMMETRY)
                && is_cs_yz(cs) && matches_gcv(gcv, {1, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::OUTFLOW && m_params.gclabel_outflow && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return m_params.outflow_label;

            if (bc == Gbc::NUMBEACH && m_params.gclabel_outflow && m_params.i10_enabled && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_yz(int cs) const
        {
            return cs == Dir::Y_POS || cs == Dir::Y_NEG
                || cs == Dir::Z_NEG || cs == Dir::Z_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_x(int cs) const
        {
            return cs == Dir::X_NEG || cs == Dir::X_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_patch_bc(int bc) const
        {
            return bc == 111 || bc == 112 || bc == 121 || bc == 122;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_gcv(int gcv, std::initializer_list<int> values) const
        {
            for (int value : values)
                if (gcv == value)
                    return true;
            return false;
        }

        Field1Params m_params{};
    };

    struct Field2BcDecision {
        struct Field2Params {
            AMREX_GPU_HOST_DEVICE
            Field2Params() noexcept = default;

            BoundaryConditionTypeLabel gclabel_v = BoundaryConditionTypeLabel::NEUMANN;
            BoundaryConditionTypeLabel orth_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel inflow_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel outflow_label = BoundaryConditionTypeLabel::NEUMANN;
            bool awa_label = false;
            bool gclabel_outflow = true;
        };

        AMREX_GPU_HOST_DEVICE
        Field2BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Field2BcDecision(const Field2Params& params)
            : m_params(params)
        {}

        AMREX_GPU_DEVICE
        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            if (gcv == 50)
                return BoundaryConditionTypeLabel::NEUMANN;

            const bool is_wavegen_numbeach_wall = (bc == Gbc::WAVEGEN || bc == Gbc::NUMBEACH || bc == Gbc::WALL);
            if (is_wavegen_numbeach_wall && is_cs_xz(cs) && matches_gcv(gcv, {11, 115}))
                return m_params.gclabel_v;

            if (is_wavegen_numbeach_wall && is_cs_xz(cs) && gcv == 111)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == Gbc::WALL && gcv == 15)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (matches_any(bc, {Gbc::NUMBEACH, Gbc::WALL}) && is_cs_y(cs) && gcv == 11)
                return m_params.orth_label;

            if (matches_any(bc, {Gbc::NUMBEACH, Gbc::WALL}) && is_cs_y(cs) && gcv == 8)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == Gbc::WAVEGEN && matches_gcv(gcv, {8, 11}))
                return m_params.inflow_label;

            if (bc == Gbc::OUTFLOW && m_params.gclabel_outflow && is_cs_xz(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::OUTFLOW && m_params.gclabel_outflow && is_cs_y(cs) && gcv == 11)
                return m_params.outflow_label;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {8, 11}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::SYMMETRY && is_cs_xz(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::SYMMETRY && is_cs_y(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::DIRICHLET_ORTH;

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_x(int cs) const
        {
            return cs == Dir::X_NEG || cs == Dir::X_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_y(int cs) const
        {
            return cs == Dir::Y_POS || cs == Dir::Y_NEG;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_xz(int cs) const
        {
            return is_cs_x(cs) || cs == Dir::Z_NEG || cs == Dir::Z_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_patch_bc(int bc) const
        {
            return bc == 111 || bc == 112 || bc == 121 || bc == 122;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_any(int value, std::initializer_list<int> list) const
        {
            for (int entry : list)
                if (value == entry)
                    return true;
            return false;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_gcv(int gcv, std::initializer_list<int> values) const
        {
            for (int value : values)
                if (gcv == value)
                    return true;
            return false;
        }

        Field2Params m_params{};
    };

    struct Field3BcDecision {
        struct Field3Params {
            AMREX_GPU_HOST_DEVICE
            Field3Params() noexcept = default;

            BoundaryConditionTypeLabel gclabel_w = BoundaryConditionTypeLabel::NEUMANN;
            BoundaryConditionTypeLabel orth_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel inflow_label = BoundaryConditionTypeLabel::DIRICHLET_ORTH;
            BoundaryConditionTypeLabel outflow_label = BoundaryConditionTypeLabel::NEUMANN;
            bool awa_label = false;
            bool gclabel_outflow = true;
            int A10 = 0;
        };

        AMREX_GPU_HOST_DEVICE
        Field3BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Field3BcDecision(const Field3Params& params)
            : m_params(params)
        {}

        AMREX_GPU_DEVICE
        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            if (gcv == 50)
                return BoundaryConditionTypeLabel::NEUMANN;

            const bool is_parallel_wall = ((bc == Gbc::NUMBEACH && !m_params.awa_label) || bc == Gbc::WALL);
            if (is_parallel_wall && is_cs_xy(cs) && matches_gcv(gcv, {12, 116}))
                return m_params.gclabel_w;

            if (is_parallel_wall && is_cs_xy(cs) && gcv == 112)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == Gbc::WALL && gcv == 16)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (is_parallel_wall && is_cs_z(cs) && gcv == 9)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (is_parallel_wall && matches_gcv(gcv, {12}) && (cs == Dir::Z_POS
                || (cs == Dir::Z_NEG && m_params.A10 == 6)))
                return m_params.orth_label;

            if (bc == Gbc::WAVEGEN && matches_gcv(gcv, {9, 12}))
                return m_params.inflow_label;

            if (bc == Gbc::OUTFLOW && m_params.gclabel_outflow && is_cs_xy(cs) && gcv == 12)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::OUTFLOW && m_params.gclabel_outflow && is_cs_z(cs) && gcv == 12)
                return m_params.outflow_label;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {9, 12}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::SYMMETRY && is_cs_xy(cs) && matches_gcv(gcv, {12, 19}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == Gbc::SYMMETRY && is_cs_z(cs) && matches_gcv(gcv, {12, 19}))
                return (m_params.A10 == 3 ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NOSLIP);

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_xy(int cs) const
        {
            return cs == Dir::X_NEG || cs == Dir::X_POS
                || cs == Dir::Y_NEG || cs == Dir::Y_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_z(int cs) const
        {
            return cs == Dir::Z_NEG || cs == Dir::Z_POS;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_patch_bc(int bc) const
        {
            return bc == 111 || bc == 112 || bc == 121 || bc == 122;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_gcv(int gcv, std::initializer_list<int> values) const
        {
            for (int value : values)
                if (gcv == value)
                    return true;
            return false;
        }

        Field3Params m_params{};
    };

    struct Field4BcDecision {
        struct Field4Params {
            AMREX_GPU_HOST_DEVICE
            Field4Params() noexcept = default;

            bool pressout_label = false;
            bool pressin_label = false;
            bool awa_label = false;
            bool gclabel_lsm_in_neumann = true;
            bool gclabel_press_in_neumann = true;
            int B77 = 0;
            int H61 = 0;
            int H62 = 0;
            int H63 = 0;
            int H64 = 0;
            int H65 = 0;
            int H66 = 0;
        };

        AMREX_GPU_HOST_DEVICE
        Field4BcDecision() = default;

        AMREX_GPU_HOST_DEVICE
        explicit Field4BcDecision(const Field4Params& params)
            : m_params(params)
        {}

        AMREX_GPU_DEVICE
        BoundaryConditionTypeLabel evaluate(int gcv, int bc, int cs) const
        {
            if (bc == 0)
                return BoundaryConditionTypeLabel::NONE;

            if (matches_gcv(gcv, {1, 30, 50, 60, 71, 72, 73, 74, 80, 150, 151, 152, 153, 154}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (gcv == 2)
            {
                if (bc != Gbc::WALL && cs != Dir::Z_NEG)
                    return BoundaryConditionTypeLabel::NEUMANN;
                return BoundaryConditionTypeLabel::NONE;
            }

            if (gcv == 40)
            {
                if ((bc == Gbc::OUTFLOW && !m_params.pressout_label)
                    || bc == Gbc::SYMMETRY
                    || (bc == Gbc::WAVEGEN && !m_params.pressin_label)
                    || (bc == Gbc::NUMBEACH && !m_params.awa_label)
                    || bc == Gbc::WALL
                    || matches_patch_bc(bc))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if ((bc == Gbc::INFLOW && !m_params.pressin_label)
                    || matches_patch_bc(bc))
                    return m_params.gclabel_press_in_neumann ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NONE;

                return BoundaryConditionTypeLabel::NONE;
            }

            if (gcv == 51 || gcv == 52 || gcv == 53 || gcv == 54)
            {
                if ((bc == Gbc::SYMMETRY || bc == Gbc::NUMBEACH || bc == Gbc::WALL)
                    || matches_patch_bc(bc))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if (((bc == Gbc::INFLOW || bc == Gbc::WAVEGEN)
                    || matches_patch_bc(bc))
                    && matches_gcv(gcv, {52, 54}))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if ((bc == Gbc::OUTFLOW)
                    || matches_patch_bc(bc))
                {
                    if (gcv == 51 || (gcv == 52 && m_params.B77 == 1) || gcv == 54)
                        return BoundaryConditionTypeLabel::NEUMANN;
                }

                if (bc == Gbc::WAVEGEN && matches_gcv(gcv, {51, 53}))
                    return m_params.gclabel_lsm_in_neumann ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NONE;

                return BoundaryConditionTypeLabel::NONE;
            }

            const bool is_symm = (bc == Gbc::SYMMETRY);
            const bool is_wavegen = (bc == Gbc::WAVEGEN);
            const bool is_numbeach = (bc == Gbc::NUMBEACH);
            const bool is_inflow = (bc == Gbc::INFLOW);
            const bool is_outflow = (bc == Gbc::OUTFLOW);
            const bool is_wall = (bc == Gbc::WALL);

            switch (gcv)
            {
                case 20:
                    if (is_symm && cs == Dir::Z_POS)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_outflow || is_symm)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 24:
                    if (is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_inflow || is_outflow || is_symm || is_wall)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 80:
                    if (heat_match(cs))
                        return BoundaryConditionTypeLabel::HEATBC;
                    return BoundaryConditionTypeLabel::NONE;
                case 49:
                    if (is_wall || is_symm)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::POTENTIAL;
                    return BoundaryConditionTypeLabel::NONE;
                case 250:
                    if ((is_numbeach || is_wall) && cs != Dir::Z_NEG)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_symm && cs != Dir::Z_POS)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 101:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::X_NEG || cs == Dir::X_POS))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::Y_POS || cs == Dir::Y_NEG || cs == Dir::Z_NEG || cs == Dir::Z_POS))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 102:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::Y_POS || cs == Dir::Y_NEG))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::X_NEG || cs == Dir::X_POS || cs == Dir::Z_NEG || cs == Dir::Z_POS))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 103:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::Z_NEG || cs == Dir::Z_POS))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == Dir::X_NEG || cs == Dir::Y_POS || cs == Dir::Y_NEG || cs == Dir::X_POS))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                default:
                    break;
            }

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_patch_bc(int bc) const
        {
            return bc == 111 || bc == 112 || bc == 121 || bc == 122;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_gbc(int bc, std::initializer_list<int> values) const
        {
            for (int value : values)
                if (bc == value)
                    return true;
            return false;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool matches_gcv(int gcv, std::initializer_list<int> values) const
        {
            for (int value : values)
                if (gcv == value)
                    return true;
            return false;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool heat_match(int cs) const
        {
            if (m_params.H61 == 1 && cs == Dir::X_NEG)
                return true;
            if (m_params.H62 == 1 && cs == Dir::Y_POS)
                return true;
            if (m_params.H63 == 1 && cs == Dir::Y_NEG)
                return true;
            if (m_params.H64 == 1 && cs == Dir::X_POS)
                return true;
            if (m_params.H65 == 1 && cs == Dir::Z_NEG)
                return true;
            if (m_params.H66 == 1 && cs == Dir::Z_POS)
                return true;
            return false;
        }

        Field4Params m_params{};
    };

    struct ConstMyExtBCFillFieldParams {
        AMREX_GPU_HOST_DEVICE
        ConstMyExtBCFillFieldParams() noexcept = default;
        ConstMyExtBCFillFieldParams(const amrex::Array<int,6>& bc_values_in,
                                const amrex::Array<amrex::Real,6>& heat_values_in,
                                bool y_dimension_exists_in, amrex_bc_func::DataLocation data_location_in)
            : bc_values(bc_values_in), heat_values(heat_values_in),
              y_dimension_exists(y_dimension_exists_in), data_location(data_location_in) {}

        //Const parameters needed for BC decision making
        const amrex::Array<int,6> bc_values = {};
        const amrex::Array<amrex::Real,6> heat_values = {};
        const bool y_dimension_exists = true;
        const DataLocation data_location = DataLocation::CELL_CENTERED;
    };
    struct MyExtBCFillFieldParams {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillFieldParams() noexcept = default;
        MyExtBCFillFieldParams(double Ui_in, double Uo_in, double dt_in, const amrex::Array<int,6>& face_labels_in)
            : Ui(Ui_in), Uo(Uo_in), dt(dt_in), face_labels(face_labels_in) {}

        double Ui = 0.0;
        double Uo = 0.0;
        double dt = 0.0;
        // Path B: pre-evaluated face labels for single-component fields.
        // Indexed as: 0=Z_NEG  1=Z_POS  2=X_NEG  3=X_POS  4=Y_NEG  5=Y_POS
        // Populated in FillDomainBoundaryImpl; avoids BCRec pointer dereference
        // in the common (numcomp==1) face path.
        amrex::Array<int,6> face_labels = {};
    };
    template <typename BCDecision>
    struct MyExtBCFillField {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillField() = delete;

        AMREX_GPU_HOST_DEVICE
        MyExtBCFillField(const ConstMyExtBCFillFieldParams const_params, const MyExtBCFillFieldParams params)
            : m_const_params(const_params), m_params(params) {}

        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                        const int dcomp, const int numcomp,
                        amrex::GeometryData const& geom, const amrex::Real time,
                        const amrex::BCRec* bcr, const int bcomp,
                        const int orig_comp) const
        {
            amrex::ignore_unused(time, orig_comp);

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
            int out_x = 0, out_y = 0, out_z = 0;

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

            if (iv[2] < dom.smallEnd(2)) {
                out_z = -1; interior[2] = dom.smallEnd(2);
            } else if (iv[2] > dom.bigEnd(2)) {
                out_z = 1; interior[2] = dom.bigEnd(2);
                if (m_const_params.data_location == DataLocation::FACE_Z) --interior[2];
            }

            if (((out_x != 0) + (out_y != 0) + (out_z != 0)) == 1)
            {
                // Pure face cell.  Encode face as an index matching face_labels:
                //   0=Z_NEG  1=Z_POS  2=X_NEG  3=X_POS  4=Y_NEG  5=Y_POS
                int face_idx, cs;
                if      (out_z == -1) { face_idx = 0; cs = Dir::Z_NEG; }
                else if (out_z ==  1) { face_idx = 1; cs = Dir::Z_POS; }
                else if (out_x == -1) { face_idx = 2; cs = Dir::X_NEG; }
                else if (out_x ==  1) { face_idx = 3; cs = Dir::X_POS; }
                else if (out_y == -1) { face_idx = 4; cs = Dir::Y_NEG; }
                else                  { face_idx = 5; cs = Dir::Y_POS; }

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
                    const int  bcr_dim = (face_idx < 2) ? 2 : (face_idx < 4) ? 0 : 1;
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

            // Edge or corner — dispatch using the already-computed out_x/out_y/out_z,
            // avoiding redundant domain-bounds re-checks (~30 comparisons eliminated).
            //
            // Behavior map (output-identical to the original if/else chain):
            //   X+Z edge (out_x!=0, out_z!=0, out_y==0)    → conditional NOSLIP/NEUMANN
            //   Y+X_MAX edge (out_y!=0, out_x>0, out_z==0) → is_corner_layer1 (preserves
            //                                                  original detection-bug behavior)
            //   true corner (all three non-zero)            → is_corner_layer1
            //   everything else (Y+X_NEG, Z+Y edges)       → NEUMANN
            const bool is_face_xz = (m_const_params.data_location == DataLocation::FACE_X
                                  || m_const_params.data_location == DataLocation::FACE_Z);
            const bool is_xz_edge = (out_x != 0 && out_z != 0 && out_y == 0);
            const bool use_corner_logic = (!is_xz_edge)
                                       && ((out_y != 0 && out_x > 0 && out_z == 0)    // Y+X_MAX (original detection bug preserved)
                                        || (out_x != 0 && out_y != 0 && out_z != 0)); // true corner

            for (int n = 0; n < numcomp; ++n)
            {
                BoundaryConditionTypeLabel label = BoundaryConditionTypeLabel::NEUMANN;

                if (is_xz_edge && is_face_xz)
                {
                    const int z_side_idx = (out_z == -1) ? Dir::Z_NEG-1 : Dir::Z_POS-1; // 4=Z_NEG, 5=Z_POS
                    if (out_x == -1) // X_NEG + Z edges
                    {
                        if ((m_const_params.bc_values[0] == Gbc::WAVEGEN
                             && static_cast<BoundaryConditionTypeLabel>(bcr[bcomp+n].lo(0)) == BoundaryConditionTypeLabel::NONE)
                         || (m_const_params.bc_values[0] == Gbc::INFLOW
                             && m_const_params.bc_values[z_side_idx] == Gbc::WALL))
                            label = BoundaryConditionTypeLabel::NOSLIP;
                    }
                    else // X_POS + Z edges
                    {
                        if (m_const_params.bc_values[1] == Gbc::NUMBEACH
                         || (m_const_params.bc_values[1] == Gbc::OUTFLOW
                             && m_const_params.bc_values[z_side_idx] == Gbc::WALL))
                            label = BoundaryConditionTypeLabel::NOSLIP;
                    }
                }
                else if (use_corner_logic)
                {
                    label = is_corner_layer1(iv, dom)
                            ? BoundaryConditionTypeLabel::NEUMANN
                            : BoundaryConditionTypeLabel::NOSLIP;
                }
                // else: NEUMANN already set above

                apply_label(label, 0, dest, iv, interior, dcomp+n, geom);
            }
        }

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
                case BoundaryConditionTypeLabel::NOSLIP:
                case BoundaryConditionTypeLabel::DIRICHLET_ORTH:
                case BoundaryConditionTypeLabel::DIRICHLET_ORTH_REFLECT:
                case BoundaryConditionTypeLabel::DIRICHLET_PARA_REFLECT:
                    dest(iv, comp) = amrex::Real(0);
                    break;
                case BoundaryConditionTypeLabel::OUTFLOWBC:
                    if(cs==Dir::X_NEG)
                        dest(iv, comp) = dest(interior, comp) - (m_params.dt/geom.CellSize(0))*m_params.Uo*(dest(interior+amrex::IntVect(1,0,0), comp)-dest(interior, comp));
                    else if(cs==Dir::X_POS)
                        dest(iv, comp) = std::max(amrex::Real(0), dest(interior, comp) - (m_params.dt/geom.CellSize(0))*m_params.Uo*(dest(interior, comp)-dest(interior+amrex::IntVect(-1,0,0), comp)));
                    else if(cs==Dir::Y_NEG)
                        dest(iv, comp) = dest(interior, comp) - (m_params.dt/geom.CellSize(1))*m_params.Uo*(dest(interior+amrex::IntVect(0,1,0), comp)-dest(interior, comp));
                    else if(cs==Dir::Y_POS)
                        dest(iv, comp) = dest(interior, comp) - (m_params.dt/geom.CellSize(1))*m_params.Uo*(dest(interior, comp)-dest(interior+amrex::IntVect(0,0,-1), comp));
                    else if(cs==Dir::Z_NEG)
                        dest(iv, comp) = dest(interior, comp) - (m_params.dt/geom.CellSize(2))*m_params.Uo*(dest(interior+amrex::IntVect(0,0,1), comp)-dest(interior, comp));
                    else if(cs==Dir::Z_POS)
                        dest(iv, comp) = dest(interior, comp) - (m_params.dt/geom.CellSize(2))*m_params.Uo*(dest(interior, comp)-dest(interior+amrex::IntVect(0,0,-1), comp));
                    break;
                case BoundaryConditionTypeLabel::POTENTIAL:
                    if(cs==Dir::X_NEG)
                        dest(iv, comp) = m_params.Ui * geom.CellSize(0) + dest(interior, comp);
                    else if(cs==Dir::X_POS)
                        dest(iv, comp) = m_params.Uo * geom.CellSize(0) + dest(interior, comp);
                    break;
                case BoundaryConditionTypeLabel::HEATBC:
                    if(cs==Dir::X_NEG)
                        dest(iv, comp) = m_const_params.heat_values[0];
                    else if(cs==Dir::X_POS)
                        dest(iv, comp) = m_const_params.heat_values[1];
                    else if(cs==Dir::Y_NEG)
                        dest(iv, comp) = m_const_params.heat_values[2];
                    else if(cs==Dir::Y_POS)
                        dest(iv, comp) = m_const_params.heat_values[3];
                    else if(cs==Dir::Z_NEG)
                        dest(iv, comp) = m_const_params.heat_values[4];
                    else if(cs==Dir::Z_POS)
                        dest(iv, comp) = m_const_params.heat_values[5];
                    break;
            }
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_corner_layer1(const amrex::IntVect& iv, const amrex::Box& box) const
        {
            int max_dist = 0;
            for(int dir=0; dir<AMREX_SPACEDIM; ++dir)
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

        const ConstMyExtBCFillFieldParams m_const_params{};
        const MyExtBCFillFieldParams m_params{};
    };
};

#endif
#endif
