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

#ifndef FIELD_AMREX_H_
#define FIELD_AMREX_H_

#include "field.h"
#include "lexer.h"
#include "amrex_bc_func.h"
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_BCUtil.H>
#include <AMReX_Gpu.H>
#include <AMReX_Vector.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <utility>
#include <vector>

class field_amrex : public field
{
public:
    virtual ~field_amrex() = default;

    /*!
     * @copydoc field_base::operator()
     */
    double& operator()(int ii, int jj, int kk, bool addOrigin = true) override;

    /*!
     * @copydoc field_base::setVal()
     */
    void setVal(double val, bool includeGhost = false) override;

    /*!
     * @copydoc field_base::FillBoundary()
     */
    void FillBoundary() override;

    void FillDomainBoundaryValue(double value, int dir, bool high) override;

protected:
    field_amrex(lexer* p, amrex_bc_func::DataLocation data_location);

    lexer *p = nullptr;
    amrex::Vector<amrex::MultiFab> mf = {};
    amrex::Vector<amrex::Vector<amrex::BCRec>> BCRecs = {};

protected:
    template <typename BCDecision>
    void FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision);

private:
    void ShiftBigBoundaryFaceInward(amrex::MultiFab& mf_in, int p_level)
    {
        int dir = -1;
        if (const_params.data_location == amrex_bc_func::DataLocation::FACE_X) dir = 0;
        else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Y) dir = 1;
        else if (const_params.data_location == amrex_bc_func::DataLocation::FACE_Z) dir = 2;

        if (dir == -1) return;

        const auto& geom = p->amrex_geometry[p_level];
        const int domain_hi = geom.Domain().bigEnd(dir);

        for (amrex::MFIter mfi(mf_in); mfi.isValid(); ++mfi)
        {
            const amrex::Box& valid_box = mfi.validbox();

            // Check if this box touches the high boundary in the specific direction
            if (valid_box.bigEnd(dir) == domain_hi)
            {
                const amrex::Box& box = mfi.fabbox();
                amrex::Array4<amrex::Real> const& arr = mf_in.array(mfi);

                int start = domain_hi;
                int end = box.bigEnd(dir) - 1;

                // Define a box collapsed to the start plane for iteration
                amrex::Box para_box = box;
                para_box.setSmall(dir, start);
                para_box.setBig(dir, start);

                if (dir == 0)
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int /*i dummy*/, int j, int k, int n)
                    {
                        for (int i = start; i <= end; ++i)
                        {
                            arr(i, j, k, n) = arr(i + 1, j, k, n);
                        }
                    });
                }
                else if (dir == 1)
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int i, int /*j dummy*/, int k, int n)
                    {
                        for (int j = start; j <= end; ++j)
                        {
                            arr(i, j, k, n) = arr(i, j + 1, k, n);
                        }
                    });
                }
                else // dir == 2
                {
                    amrex::ParallelFor(para_box, mf_in.nComp(),
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k dummy*/, int n)
                    {
                        for (int k = start; k <= end; ++k)
                        {
                            arr(i, j, k, n) = arr(i, j, k + 1, n);
                        }
                    });
                }
            }
        }
    }

    const amrex_bc_func::ConstMyExtBCFillFieldParams const_params = {};
};

namespace field_amrex_detail
{
    using Label = amrex_bc_func::BoundaryConditionTypeLabel;

    enum class OutflowAxis { U, V, W };

    inline Label base_parallel_label(int B20)
    {
        switch (B20)
        {
            case 1: return Label::NEUMANN;
            case 2: case 4: return Label::NOSLIP;
            case 3: return Label::DIRICHLET_ORTH;
            default: return Label::NEUMANN;
        }
    }

    inline std::pair<Label, Label> compute_parallel_u_v(lexer* p)
    {
        Label u = base_parallel_label(p->B20);
        Label v = base_parallel_label(p->B20);
        if (p->B23 == 2)
        {
            u = Label::DIRICHLET_PARA_REFLECT;
            v = Label::DIRICHLET_PARA_REFLECT;
        }

        if (p->A217 == 1 && p->A10 == 2)
        {
            u = Label::NEUMANN;
            v = Label::NEUMANN;
        }
        else if (p->A217 == 2 && p->A10 == 2)
        {
            u = Label::NOSLIP;
            v = Label::NOSLIP;
        }

        return {u, v};
    }

    inline Label compute_parallel_w(lexer* p)
    {
        if (p->B23 == 2)
            return Label::DIRICHLET_PARA_REFLECT;
        return base_parallel_label(p->B20);
    }

    inline Label compute_orth_label(lexer* p)
    {
        return (p->B23 == 2) ? Label::DIRICHLET_ORTH_REFLECT : Label::DIRICHLET_ORTH;
    }

    inline Label compute_inflow_label(lexer* p)
    {
        if (p->I230 > 0 || p->B98 >= 3 || p->B60 > 0)
            return Label::NONE;
        return Label::DIRICHLET_ORTH;
    }

    inline Label compute_outflow_label(lexer* p, OutflowAxis axis)
    {
        switch (p->B75)
        {
            case 1:
                return Label::NEUMANN;
            case 2:
                return Label::OUTFLOWBC;
            case 3:
                return (axis == OutflowAxis::U) ? Label::NONE : Label::OUTFLOWBC;
            default:
                return Label::NEUMANN;
        }
    }

    inline bool compute_awa_label(lexer* p)
    {
        return p->B99 >= 3;
    }

    inline bool compute_gclabel_outflow(lexer* p)
    {
        return !(p->B60 == 3 || p->B60 == 4);
    }

    inline bool compute_i10_enabled(lexer* p)
    {
        return p->I10 == 1;
    }
}

template <typename BCDecision>
void field_amrex::FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision)
{
    amrex_bc_func::MyExtBCFillFieldParams params;
    params.Ui = p->Ui;
    params.Uo = p->Uo;
    params.dt = p->dt;

    // Update BCRecs (and local references)
    // Note: Assuming nComp is handled by resizing BCRecs correctly in constructor
    const int x_neg = 1; // Dir::X_NEG (1) in amrex_bc_func::Dir
    const int x_pos = 4; // Dir::X_POS (4)
    const int y_neg = 3; // Dir::Y_NEG (3)
    const int y_pos = 2; // Dir::Y_POS (2)
    const int z_neg = 5; // Dir::Z_NEG (5)
    const int z_pos = 6; // Dir::Z_POS (6)
    // Face indices in bc_values: 0=X-, 1=X+, 2=Y-, 3=Y+, 4=Z-, 5=Z+

    LevelLOOP
    {
        for (int n = 0; n < mf[p->level].nComp(); ++n)
        {
            if (BCRecs[p->level].size() <= n) continue;

            auto& bc = BCRecs[p->level][n];

            int bc_code_1 = const_params.bc_values[0];
            auto label_1 = bc_decision.evaluate(gcv, bc_code_1, x_neg);
            bc.setLo(0, static_cast<int>(label_1));

            int bc_code_2 = const_params.bc_values[1];
            auto label_2 = bc_decision.evaluate(gcv, bc_code_2, x_pos);
            bc.setHi(0, static_cast<int>(label_2));

            int bc_code_3 = const_params.bc_values[2];
            auto label_3 = bc_decision.evaluate(gcv, bc_code_3, y_neg);
            bc.setLo(1, static_cast<int>(label_3));

            int bc_code_4 = const_params.bc_values[3];
            auto label_4 = bc_decision.evaluate(gcv, bc_code_4, y_pos);
            bc.setHi(1, static_cast<int>(label_4));

            int bc_code_5 = const_params.bc_values[4];
            auto label_5 = bc_decision.evaluate(gcv, bc_code_5, z_neg);
            bc.setLo(2, static_cast<int>(label_5));

            int bc_code_6 = const_params.bc_values[5];
            auto label_6 = bc_decision.evaluate(gcv, bc_code_6, z_pos);
            bc.setHi(2, static_cast<int>(label_6));
        }
    }

    LevelLOOP
    {
        if(p->level==0)
        {
            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> bf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> physbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], bf);

            amrex::FillPatchSingleLevel(mf[p->level], amrex::Real(p->simtime),
                                        {&(mf[p->level])}, {amrex::Real(p->simtime)},
                                        0, 0, mf[p->level].nComp(), p->amrex_geometry[p->level], physbcf, 0);
        }
        else
        {
            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> cbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> cphysbcf(
                p->amrex_geometry[p->level-1], BCRecs[p->level-1], cbf);

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> fbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{const_params, params});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> fphysbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], fbf);

            amrex::CellConservativeLinear mapper;
            const amrex::IntVect ratio = p->ref_ratio * amrex::IntVect::TheUnitVector();

            amrex::FillPatchTwoLevels(mf[p->level], amrex::Real(p->simtime),
                                        {&(mf[p->level-1])}, {amrex::Real(p->simtime)},
                                        {&(mf[p->level])}, {amrex::Real(p->simtime)},
                                        0, 0, mf[p->level].nComp(), p->amrex_geometry[p->level-1], p->amrex_geometry[p->level],
                                        cphysbcf, 0,
                                        fphysbcf, 0, // second one?
                                        ratio, // refinement ratio
                                        &mapper, // spatial interpolater
                                        BCRecs[p->level], 0);
        }

        ShiftBigBoundaryFaceInward(mf[p->level], p->level);
    }
}

#endif
