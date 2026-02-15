/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef AMREX_BC_FUNC_H_
#define AMREX_BC_FUNC_H_

#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <initializer_list>
#include <type_traits>

class amrex_bc_func
{
public:
    amrex_bc_func() = default;
    virtual ~amrex_bc_func() = default;

    enum class BoundaryConditionTypeLabel : int { NONE = 0, DIRICHLET_ORTH = 1, NEUMANN = 4, NOSLIP = 5, OUTFLOWBC = 6, SOMMERFELD = 7,
                                        POTENTIAL = 8, DIRICHLET_ORTH_REFLECT = 11, DIRICHLET_PARA_REFLECT = 12,
                                        NEUMANN_X = 14, NEUMANN_HX = 41, NEUMANN_HY = 42, HEATBC = 61 };
private:
    enum class Gbc : int { INFLOW = 1, OUTFLOW = 2, SYMMETRY = 3, WAVEGEN = 6, NUMBEACH = 7, WALL = 21 };
    enum class Dir : int { X_NEG = 1, X_POS = 4, Y_NEG = 3, Y_POS = 2, Z_NEG = 5, Z_POS = 6 };
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

            const bool is_parallel_wall = ((bc == static_cast<int>(Gbc::NUMBEACH) && !m_params.awa_label) || bc == static_cast<int>(Gbc::WALL));
            if (is_parallel_wall && is_cs_yz(cs) && matches_gcv(gcv, {1, 10, 114}))
                return m_params.gclabel_u;

            if (is_parallel_wall && is_cs_yz(cs) && gcv == 110)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == static_cast<int>(Gbc::WALL) && gcv == 14)
                return BoundaryConditionTypeLabel::NEUMANN;

            if ((bc == static_cast<int>(Gbc::SYMMETRY) || bc == static_cast<int>(Gbc::WALL)) && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return m_params.orth_label;

            if (bc == static_cast<int>(Gbc::WAVEGEN) && is_cs_x(cs) && matches_gcv(gcv, {1, 7, 10}))
                return m_params.inflow_label;

            if (bc == static_cast<int>(Gbc::WALL) && is_cs_x(cs) && gcv == 7)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {1, 7, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (((bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow) || bc == static_cast<int>(Gbc::SYMMETRY))
                && is_cs_yz(cs) && matches_gcv(gcv, {1, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return m_params.outflow_label;

            if (bc == static_cast<int>(Gbc::NUMBEACH) && m_params.gclabel_outflow && m_params.i10_enabled && is_cs_x(cs) && matches_gcv(gcv, {1, 10}))
                return BoundaryConditionTypeLabel::NEUMANN;

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_yz(int cs) const
        {
            return cs == static_cast<int>(Dir::Y_POS) || cs == static_cast<int>(Dir::Y_NEG)
                || cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_x(int cs) const
        {
            return cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::X_POS);
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

            const bool is_wavegen_numbeach_wall = (bc == static_cast<int>(Gbc::WAVEGEN) || bc == static_cast<int>(Gbc::NUMBEACH) || bc == static_cast<int>(Gbc::WALL));
            if (is_wavegen_numbeach_wall && is_cs_xz(cs) && matches_gcv(gcv, {11, 115}))
                return m_params.gclabel_v;

            if (is_wavegen_numbeach_wall && is_cs_xz(cs) && gcv == 111)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == static_cast<int>(Gbc::WALL) && gcv == 15)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (matches_any(bc, {static_cast<int>(Gbc::NUMBEACH), static_cast<int>(Gbc::WALL)}) && is_cs_y(cs) && gcv == 11)
                return m_params.orth_label;

            if (matches_any(bc, {static_cast<int>(Gbc::NUMBEACH), static_cast<int>(Gbc::WALL)}) && is_cs_y(cs) && gcv == 8)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == static_cast<int>(Gbc::WAVEGEN) && matches_gcv(gcv, {8, 11}))
                return m_params.inflow_label;

            if (bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow && is_cs_xz(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow && is_cs_y(cs) && gcv == 11)
                return m_params.outflow_label;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {8, 11}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::SYMMETRY) && is_cs_xz(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::SYMMETRY) && is_cs_y(cs) && gcv == 11)
                return BoundaryConditionTypeLabel::DIRICHLET_ORTH;

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_x(int cs) const
        {
            return cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::X_POS);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_y(int cs) const
        {
            return cs == static_cast<int>(Dir::Y_POS) || cs == static_cast<int>(Dir::Y_NEG);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_xz(int cs) const
        {
            return is_cs_x(cs) || cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS);
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

            const bool is_parallel_wall = ((bc == static_cast<int>(Gbc::NUMBEACH) && !m_params.awa_label) || bc == static_cast<int>(Gbc::WALL));
            if (is_parallel_wall && is_cs_xy(cs) && matches_gcv(gcv, {12, 116}))
                return m_params.gclabel_w;

            if (is_parallel_wall && is_cs_xy(cs) && gcv == 112)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (bc == static_cast<int>(Gbc::WALL) && gcv == 16)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (is_parallel_wall && is_cs_z(cs) && gcv == 9)
                return BoundaryConditionTypeLabel::NOSLIP;

            if (is_parallel_wall && matches_gcv(gcv, {12}) && (cs == static_cast<int>(Dir::Z_POS)
                || (cs == static_cast<int>(Dir::Z_NEG) && m_params.A10 == 6)))
                return m_params.orth_label;

            if (bc == static_cast<int>(Gbc::WAVEGEN) && matches_gcv(gcv, {9, 12}))
                return m_params.inflow_label;

            if (bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow && is_cs_xy(cs) && gcv == 12)
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::OUTFLOW) && m_params.gclabel_outflow && is_cs_z(cs) && gcv == 12)
                return m_params.outflow_label;

            if (matches_patch_bc(bc) && matches_gcv(gcv, {9, 12}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::SYMMETRY) && is_cs_xy(cs) && matches_gcv(gcv, {12, 19}))
                return BoundaryConditionTypeLabel::NEUMANN;

            if (bc == static_cast<int>(Gbc::SYMMETRY) && is_cs_z(cs) && matches_gcv(gcv, {12, 19}))
                return (m_params.A10 == 3 ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NOSLIP);

            return BoundaryConditionTypeLabel::NONE;
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_xy(int cs) const
        {
            return cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::X_POS)
                || cs == static_cast<int>(Dir::Y_NEG) || cs == static_cast<int>(Dir::Y_POS);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_cs_z(int cs) const
        {
            return cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS);
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
                if (bc != 5 && bc != static_cast<int>(Gbc::WALL) && cs != static_cast<int>(Dir::Z_NEG))
                    return BoundaryConditionTypeLabel::NEUMANN;
                return BoundaryConditionTypeLabel::NONE;
            }

            if (gcv == 40)
            {
                if ((bc == static_cast<int>(Gbc::OUTFLOW) && !m_params.pressout_label)
                    || bc == static_cast<int>(Gbc::SYMMETRY)
                    || (bc == static_cast<int>(Gbc::WAVEGEN) && !m_params.pressin_label)
                    || (bc == static_cast<int>(Gbc::NUMBEACH) && !m_params.awa_label)
                    || bc == static_cast<int>(Gbc::WALL)
                    || matches_patch_bc(bc))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if ((bc == static_cast<int>(Gbc::INFLOW) && !m_params.pressin_label)
                    || matches_patch_bc(bc))
                    return m_params.gclabel_press_in_neumann ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NONE;

                return BoundaryConditionTypeLabel::NONE;
            }

            if (gcv == 51 || gcv == 52 || gcv == 53 || gcv == 54)
            {
                if ((bc == static_cast<int>(Gbc::SYMMETRY) || bc == static_cast<int>(Gbc::NUMBEACH) || bc == static_cast<int>(Gbc::WALL))
                    || matches_patch_bc(bc))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if (((bc == static_cast<int>(Gbc::INFLOW) || bc == static_cast<int>(Gbc::WAVEGEN))
                    || matches_patch_bc(bc))
                    && matches_gcv(gcv, {52, 54}))
                    return BoundaryConditionTypeLabel::NEUMANN;

                if ((bc == static_cast<int>(Gbc::OUTFLOW))
                    || matches_patch_bc(bc))
                {
                    if (gcv == 51 || (gcv == 52 && m_params.B77 == 1) || gcv == 54)
                        return BoundaryConditionTypeLabel::NEUMANN;
                }

                if (bc == static_cast<int>(Gbc::WAVEGEN) && matches_gcv(gcv, {51, 53}))
                    return m_params.gclabel_lsm_in_neumann ? BoundaryConditionTypeLabel::NEUMANN : BoundaryConditionTypeLabel::NONE;

                return BoundaryConditionTypeLabel::NONE;
            }

            const bool is_symm = (bc == static_cast<int>(Gbc::SYMMETRY));
            const bool is_wavegen = (bc == static_cast<int>(Gbc::WAVEGEN));
            const bool is_numbeach = (bc == static_cast<int>(Gbc::NUMBEACH));
            const bool is_inflow = (bc == static_cast<int>(Gbc::INFLOW));
            const bool is_outflow = (bc == static_cast<int>(Gbc::OUTFLOW));
            const bool is_wall = (bc == static_cast<int>(Gbc::WALL));

            switch (gcv)
            {
                case 20:
                    if (is_symm && cs == static_cast<int>(Dir::Z_POS))
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
                    if ((is_numbeach || is_wall) && cs != static_cast<int>(Dir::Z_NEG))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_symm && cs != static_cast<int>(Dir::Z_POS))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 101:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::X_POS)))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::Y_POS) || cs == static_cast<int>(Dir::Y_NEG) || cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS)))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 102:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::Y_POS) || cs == static_cast<int>(Dir::Y_NEG)))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::X_POS) || cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS)))
                        return BoundaryConditionTypeLabel::NEUMANN;
                    if (is_inflow || is_outflow || is_wavegen || is_numbeach)
                        return BoundaryConditionTypeLabel::NEUMANN;
                    return BoundaryConditionTypeLabel::NONE;
                case 103:
                    if (is_wall)
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::Z_NEG) || cs == static_cast<int>(Dir::Z_POS)))
                        return BoundaryConditionTypeLabel::NOSLIP;
                    if (is_symm && (cs == static_cast<int>(Dir::X_NEG) || cs == static_cast<int>(Dir::Y_POS) || cs == static_cast<int>(Dir::Y_NEG) || cs == static_cast<int>(Dir::X_POS)))
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
            if (m_params.H61 == 1 && cs == static_cast<int>(Dir::X_NEG))
                return true;
            if (m_params.H62 == 1 && cs == static_cast<int>(Dir::Y_POS))
                return true;
            if (m_params.H63 == 1 && cs == static_cast<int>(Dir::Y_NEG))
                return true;
            if (m_params.H64 == 1 && cs == static_cast<int>(Dir::X_POS))
                return true;
            if (m_params.H65 == 1 && cs == static_cast<int>(Dir::Z_NEG))
                return true;
            if (m_params.H66 == 1 && cs == static_cast<int>(Dir::Z_POS))
                return true;
            return false;
        }

        Field4Params m_params{};
    };

    struct ConstMyExtBCFillFieldParams {
        AMREX_GPU_HOST_DEVICE
        ConstMyExtBCFillFieldParams() noexcept = default;

        //Const parameters needed for BC decision making
        const amrex::Array<int,6> bc_values = {};
        const amrex::Array<amrex::Real,6> heat_values = {};
        const int margin = 0;
        const int orderdir = 3;
        const bool y_dimension_exists = true;
        const double gamma = 0;
    };
    struct MyExtBCFillFieldParams {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillFieldParams() noexcept = default;

        //Const parameters needed for BC decision making
        double Ui = 0.0;
        double Uo = 0.0;
        double dt = 0.0;
    };
    template <typename BCDecision>
    struct MyExtBCFillField {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillField() = delete;

        AMREX_GPU_HOST_DEVICE
        MyExtBCFillField(const ConstMyExtBCFillFieldParams const_params, const MyExtBCFillFieldParams params,
                const int gcv, const BCDecision& decision,
                const amrex::Box* boxes, const int num_boxes)
            : m_const_params(const_params), m_params(params), m_gcv(gcv), m_bc_decision(decision), m_boxes(boxes), m_num_boxes(num_boxes) {}

        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                        const int dcomp, const int numcomp,
                        amrex::GeometryData const& geom, const amrex::Real time,
                        const amrex::BCRec* bcr, const int bcomp,
                        const int orig_comp) const
        {
            amrex::ignore_unused(time, bcr, bcomp, orig_comp);

            if (!m_const_params.y_dimension_exists && iv[1] != 0)
                return;

            BoundaryConditionTypeLabel label = BoundaryConditionTypeLabel::NONE;
            const amrex::Box* matched_box = nullptr;
            int face_for_bc = 0;
            bool is_corner = false;
            int cs = 0;

            for (int idx = 0; idx < m_num_boxes; ++idx)
            {
                const amrex::Box& box = m_boxes[idx];
                int face = detect_face(iv, box);
                if (face > 0)
                {
                    int bc_code = m_const_params.bc_values[face-1];
                    if (bc_code == 0)
                        continue;
                    cs = cs_from_face(face);
                    label = m_bc_decision.evaluate(m_gcv, bc_code, cs);
                    if (label != BoundaryConditionTypeLabel::NONE)
                    {
                        matched_box = &box;
                        face_for_bc = face;
                        break;
                    }
                }
                else if (is_corner_layer1(iv, box))
                {
                    // label = evaluate_corner(iv, box, face_for_bc);
                    // if (label != BoundaryConditionTypeLabel::NONE)
                    // {
                        matched_box = &box;
                        is_corner = true;
                        label = BoundaryConditionTypeLabel::NEUMANN;
                        break;
                    // }
                }
            }

            if (label == BoundaryConditionTypeLabel::NONE || matched_box == nullptr)
                return;

            const amrex::Box dom = geom.Domain();
            amrex::IntVect interior = iv;
            for(int dir=0; dir<AMREX_SPACEDIM; ++dir)
            {
                if(interior[dir] < dom.smallEnd(dir))
                    interior[dir] = dom.smallEnd(dir);
                else if(interior[dir] > dom.bigEnd(dir))
                    interior[dir] = dom.bigEnd(dir);
            }

            if (!is_corner && !is_within_margin(iv, *matched_box, face_for_bc))
                return;

            for(int n=0; n<numcomp; ++n)
            {
                switch (label)
                {
                    case BoundaryConditionTypeLabel::DIRICHLET_ORTH:
                        // ToDo
                        break;
                    case BoundaryConditionTypeLabel::NEUMANN:
                    default:
                        dest(iv, dcomp+n) = dest(interior, dcomp+n);
                        break;
                    case BoundaryConditionTypeLabel::NOSLIP:
                        dest(iv, dcomp+n) = amrex::Real(0);
                        break;
                    case BoundaryConditionTypeLabel::OUTFLOWBC:
                        if(cs==static_cast<int>(Dir::X_NEG))
                            dest(iv, dcomp+n) = dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(0))*m_params.Uo*(dest(interior+amrex::IntVect(1,0,0), dcomp+n)-dest(interior, dcomp+n));
                        else if(cs==static_cast<int>(Dir::X_POS))
                            dest(iv, dcomp+n) = MAX(amrex::Real(0), dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(0))*m_params.Uo*(dest(interior, dcomp+n)-dest(interior+amrex::IntVect(-1,0,0), dcomp+n)));
                        else if(cs==static_cast<int>(Dir::Y_NEG))
                            dest(iv, dcomp+n) = dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(1))*m_params.Uo*(dest(interior+amrex::IntVect(0,1,0), dcomp+n)-dest(interior, dcomp+n));
                        else if(cs==static_cast<int>(Dir::Y_POS))
                            dest(iv, dcomp+n) = dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(1))*m_params.Uo*(dest(interior, dcomp+n)-dest(interior+amrex::IntVect(0,0,-1), dcomp+n));
                        else if(cs==static_cast<int>(Dir::Z_NEG))
                            dest(iv, dcomp+n) = dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(2))*m_params.Uo*(dest(interior+amrex::IntVect(0,0,1), dcomp+n)-dest(interior, dcomp+n));
                        else if(cs==static_cast<int>(Dir::Z_POS))
                            dest(iv, dcomp+n) = dest(interior, dcomp+n) - (m_params.dt/geom.CellSize(2))*m_params.Uo*(dest(interior, dcomp+n)-dest(interior+amrex::IntVect(0,0,-1), dcomp+n));
                        break;
                    case BoundaryConditionTypeLabel::POTENTIAL:
                        if(cs==static_cast<int>(Dir::X_NEG))
                            dest(iv, dcomp+n) = m_params.Ui * geom.CellSize(0) + dest(interior, dcomp+n);
                        else if(cs==static_cast<int>(Dir::X_POS))
                            dest(iv, dcomp+n) = m_params.Uo * geom.CellSize(0) + dest(interior, dcomp+n);
                        break;
                    case BoundaryConditionTypeLabel::DIRICHLET_ORTH_REFLECT:
                        // ToDo
                        break;
                    case BoundaryConditionTypeLabel::DIRICHLET_PARA_REFLECT:
                        // ToDo
                        break;
                    case BoundaryConditionTypeLabel::HEATBC:
                        if(cs==static_cast<int>(Dir::X_NEG))
                            dest(iv, dcomp+n) = m_const_params.heat_values[0];
                        else if(cs==static_cast<int>(Dir::X_POS))
                            dest(iv, dcomp+n) = m_const_params.heat_values[1];
                        else if(cs==static_cast<int>(Dir::Y_NEG))
                            dest(iv, dcomp+n) = m_const_params.heat_values[2];
                        else if(cs==static_cast<int>(Dir::Y_POS))
                            dest(iv, dcomp+n) = m_const_params.heat_values[3];
                        else if(cs==static_cast<int>(Dir::Z_NEG))
                            dest(iv, dcomp+n) = m_const_params.heat_values[4];
                        else if(cs==static_cast<int>(Dir::Z_POS))
                            dest(iv, dcomp+n) = m_const_params.heat_values[5];
                        break;
                }
            }
        }

    private:
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int detect_face(const amrex::IntVect& iv, const amrex::Box& dom) const
        {
            if(iv[0] < dom.smallEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                return 1;
            else if(iv[0] > dom.bigEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                return 2;
            else if(iv[1] < dom.smallEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                return 3;
            else if(iv[1] > dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                return 4;
            else if(iv[2] < dom.smallEnd(2) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                return 5;
            else if(iv[2] > dom.bigEnd(2) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                return 6;
            else
                return 0;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        BoundaryConditionTypeLabel evaluate_corner(const amrex::IntVect& iv, const amrex::Box& box, int& face_out) const
        {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                if (iv[dir] < box.smallEnd(dir))
                {
                    int face = face_for_dir(dir, false);
                    face_out = face;
                    int bc_code = m_const_params.bc_values[face-1];
                    if (bc_code == 0) continue;
                    return m_bc_decision.evaluate(m_gcv, bc_code, cs_from_face(face));
                }
                else if (iv[dir] > box.bigEnd(dir))
                {
                    int face = face_for_dir(dir, true);
                    face_out = face;
                    int bc_code = m_const_params.bc_values[face-1];
                    if (bc_code == 0) continue;
                    return m_bc_decision.evaluate(m_gcv, bc_code, cs_from_face(face));
                }
            }
            return BoundaryConditionTypeLabel::NONE;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_corner_layer1(const amrex::IntVect& iv, const amrex::Box& box) const
        {
            int outside_dims = 0;
            int max_dist = 0;
            for(int dir=0; dir<AMREX_SPACEDIM; ++dir) {
                int d = 0;
                if (iv[dir] < box.smallEnd(dir)) {
                    d = box.smallEnd(dir) - iv[dir];
                    outside_dims++;
                }
                else if (iv[dir] > box.bigEnd(dir)) {
                    d = iv[dir] - box.bigEnd(dir);
                    outside_dims++;
                }
                if (d > max_dist) max_dist = d;
            }

            return (outside_dims >= 2) && (max_dist == 1);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool is_within_margin(const amrex::IntVect& iv, const amrex::Box& box, int face) const
        {
            if (face <= 0 || m_const_params.margin <= 0)
                return false;

            int dir = (face < 3) ? 0 : (face < 5 ? 1 : 2);
            bool high = (face % 2) == 0;
            int boundary = high ? box.bigEnd(dir) : box.smallEnd(dir);
            int dist = high ? iv[dir] - boundary : boundary - iv[dir];
            return (dist >= 1) && (dist <= m_const_params.margin);
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int cs_from_face(int face) const
        {
            constexpr int map[] = {0, 1, 4, 3, 2, 5, 6};
            return map[face];
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int face_for_dir(int dir, bool high) const
        {
            if (dir == 0)
                return high ? 2 : 1;
            else if (dir == 1)
                return high ? 4 : 3;
            else
                return high ? 6 : 5;
        }

        const ConstMyExtBCFillFieldParams m_const_params{};
        const MyExtBCFillFieldParams m_params{};
        BCDecision m_bc_decision;
        const amrex::Box* m_boxes = nullptr;
        const int m_num_boxes = 0;
        const int m_gcv = 0;
    };
};

#endif