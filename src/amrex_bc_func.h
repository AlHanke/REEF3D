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
    enum class DataLocation : unsigned int { CELL_CENTERED = 0, FACE_X = 1, FACE_Y = 2, FACE_Z = 3 };
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
        ConstMyExtBCFillFieldParams(const amrex::Array<int,6>& bc_values_in,
                                const amrex::Array<amrex::Real,6>& heat_values_in,
                                bool y_dimension_exists_in, unsigned int data_location_in)
            : bc_values(bc_values_in), heat_values(heat_values_in),
              y_dimension_exists(y_dimension_exists_in), data_location(data_location_in) {}

        //Const parameters needed for BC decision making
        const amrex::Array<int,6> bc_values = {};
        const amrex::Array<amrex::Real,6> heat_values = {};
        const bool y_dimension_exists = true;
        const unsigned int data_location = static_cast<unsigned int>(DataLocation::CELL_CENTERED);
    };
    struct MyExtBCFillFieldParams {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFillFieldParams() noexcept = default;
        MyExtBCFillFieldParams(double Ui_in, double Uo_in, double dt_in)
            : Ui(Ui_in), Uo(Uo_in), dt(dt_in) {}

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

            // Determine primary face/direction
            // Note: This logic assumes simple orthogonal exterior check.
            // Only one direction is picked for BC application if in corner?
            // AMReX typically invokes for specific fill regions.
            // MyExtBCFillField needs to decide which BC applies.

            // Replicate detect_face logic but relative to Domain
            int face = 0;
            int edge = 0;
            if(iv[2] < dom.smallEnd(2) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                face = 5; // Z Min
            else if(iv[2] > dom.bigEnd(2) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                face = 6; // Z Max
            else if(iv[0] < dom.smallEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                face = 1; // X Min
            else if(iv[0] > dom.bigEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                face = 2; // X Max
            else if(iv[1] < dom.smallEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                face = 3; // Y Min
            else if(iv[1] > dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                face = 4; // Y Max
            else if(iv[0] < dom.smallEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] < dom.smallEnd(2))
                edge = 1; // X Min, Z Min
            else if(iv[0] < dom.smallEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] > dom.bigEnd(2))
                edge = 2; // X Min, Z Max
            else if(iv[0] > dom.bigEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] < dom.smallEnd(2))
                edge = 3; // X Max, Z Min
            else if(iv[0] > dom.bigEnd(0) && iv[1] >= dom.smallEnd(1) && iv[1] <= dom.bigEnd(1) && iv[2] > dom.bigEnd(2))
                edge = 4; // X Max, Z Max
            else if(iv[1] < dom.smallEnd(1) && iv[0] < dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                edge = 5; // Y Min, X Min
            else if(iv[1] < dom.smallEnd(1) && iv[0] > dom.bigEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                edge = 6; // Y Min, X Max
            else if(iv[1] > dom.bigEnd(1) && iv[0] < dom.smallEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                edge = 7; // Y Max, X Min
            else if(iv[1] >  dom.bigEnd(1) && iv[0] > dom.bigEnd(0) && iv[0] <= dom.bigEnd(0) && iv[2] >= dom.smallEnd(2) && iv[2] <= dom.bigEnd(2))
                edge = 8; // Y Max, X Max
            else if(iv[2] < dom.smallEnd(2) && iv[1] < dom.smallEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                edge = 9; // Z Min, Y Min
            else if(iv[2] < dom.smallEnd(2) && iv[1] > dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                edge = 10; // Z Min, Y Max
            else if(iv[2] > dom.bigEnd(2) && iv[1] < dom.smallEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                edge = 11; // Z Max, Y Min
            else if(iv[2] > dom.bigEnd(2) && iv[1] > dom.bigEnd(1) && iv[0] >= dom.smallEnd(0) && iv[0] <= dom.bigEnd(0))
                edge = 12; // Z Max, Y Max

            // Compute interior coordinate
            amrex::IntVect interior = iv;
            for(int dir=0; dir<AMREX_SPACEDIM; ++dir)
            {
                if(interior[dir] < dom.smallEnd(dir))
                    interior[dir] = dom.smallEnd(dir);
                else if(interior[dir] > dom.bigEnd(dir))
                {
                    interior[dir] = dom.bigEnd(dir);
                    if(m_const_params.data_location == dir +1)
                        interior[dir] -= 1;
                }
            }

            for(int n=0; n<numcomp; ++n)
            {
                BoundaryConditionTypeLabel label = BoundaryConditionTypeLabel::NONE;
                int cs = 0;

                // Map face to Dirk/CS and get label from BCRec
                // Face mapping: 1->X_NEG(1), 2->X_POS(4), 3->Y_NEG(3), 4->Y_POS(2), 5->Z_NEG(5), 6->Z_POS(6)
                const int bcrec_idx = bcomp + n;
                switch(face) {
                    case 1:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(0));
                        cs = static_cast<int>(Dir::X_NEG);
                        break;
                    case 2:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(0));
                        cs = static_cast<int>(Dir::X_POS);
                        break;
                    case 3:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(1));
                        cs = static_cast<int>(Dir::Y_NEG);
                        break;
                    case 4:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(1));
                        cs = static_cast<int>(Dir::Y_POS);
                        break;
                    case 5:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(2));
                        cs = static_cast<int>(Dir::Z_NEG);
                        break;
                    case 6:
                        label = static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(2));
                        cs = static_cast<int>(Dir::Z_POS);
                        break;
                    case 0:
                    switch (edge)
                    {
                        case 1:
                            if(bcr[bcrec_idx].lo(0) == bcr[bcrec_idx].lo(2) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(0)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 2:
                            if(bcr[bcrec_idx].lo(0) == bcr[bcrec_idx].hi(2) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(0)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 3:
                            if(bcr[bcrec_idx].hi(0) == bcr[bcrec_idx].lo(2) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(0)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 4:
                            if(bcr[bcrec_idx].hi(0) == bcr[bcrec_idx].hi(2) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(0)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 5:
                            if(bcr[bcrec_idx].lo(1) == bcr[bcrec_idx].lo(0) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(1)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 6:
                            if(bcr[bcrec_idx].lo(1) == bcr[bcrec_idx].hi(0) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(1)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 7:
                            if(bcr[bcrec_idx].hi(1) == bcr[bcrec_idx].lo(0) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(1)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 8:
                            if(bcr[bcrec_idx].hi(1) == bcr[bcrec_idx].hi(0) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(1)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 9:
                            if(bcr[bcrec_idx].lo(2) == bcr[bcrec_idx].lo(1) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(2)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 10:
                            if(bcr[bcrec_idx].lo(2) == bcr[bcrec_idx].hi(1) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].lo(2)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 11:
                            if(bcr[bcrec_idx].hi(2) == bcr[bcrec_idx].lo(1) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(2)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 12:
                            if(bcr[bcrec_idx].hi(2) == bcr[bcrec_idx].hi(1) && static_cast<BoundaryConditionTypeLabel>(bcr[bcrec_idx].hi(2)) == BoundaryConditionTypeLabel::NEUMANN)
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            else
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            break;
                        case 0:
                            if(is_corner_layer1(iv, dom))
                                label = BoundaryConditionTypeLabel::NOSLIP;
                            else
                                label = BoundaryConditionTypeLabel::NEUMANN;
                            break;
                    }
                }

                switch (label)
                {
                    case BoundaryConditionTypeLabel::NONE:
                        break;
                    case BoundaryConditionTypeLabel::NEUMANN:
                    default:
                        dest(iv, dcomp+n) = dest(interior, dcomp+n);
                        break;
                    case BoundaryConditionTypeLabel::NOSLIP:
                    case BoundaryConditionTypeLabel::DIRICHLET_ORTH:
                    case BoundaryConditionTypeLabel::DIRICHLET_ORTH_REFLECT:
                    case BoundaryConditionTypeLabel::DIRICHLET_PARA_REFLECT:
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
