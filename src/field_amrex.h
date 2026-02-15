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

    double& operator()(int ii, int jj, int kk) override;

    void setVal(double val, bool includeGhost = false) override;

    void FillBoundary() override;

protected:
    field_amrex(lexer* p);

    lexer *p;
    amrex::Vector<amrex::MultiFab> mf;
    amrex::Vector<amrex::Vector<amrex::BCRec>> BCRecs;
private:
    const amrex::Array<int,2*AMREX_SPACEDIM>face_bc_values;
    const amrex::Array<int,2*AMREX_SPACEDIM>face_heat_values;

protected:
    template <typename BCDecision>
    void FillDomainBoundaryImpl(int gcv, const BCDecision& bc_decision);
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
    LevelLOOP
        if(p->level==0)
        {
            amrex::Vector<amrex::Box> loc_boxes;
            for (amrex::MFIter mfi(mf[p->level]); mfi.isValid(); ++mfi)
            {
                loc_boxes.push_back(mfi.validbox());
            }

            amrex::Gpu::DeviceVector<amrex::Box> device_boxes(loc_boxes.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, loc_boxes.begin(), loc_boxes.end(), device_boxes.begin());

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> bf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{face_bc_values, face_heat_values, gcv, p->margin, static_cast<bool>(p->j_dir), bc_decision,
                                                        device_boxes.data(), static_cast<int>(device_boxes.size())});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> physbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], bf);

            amrex::FillPatchSingleLevel(mf[p->level], amrex::Real(p->simtime),
                                        {&(mf[p->level])}, {amrex::Real(p->simtime)},
                                        0, 0, mf[p->level].nComp(), p->amrex_geometry[p->level], physbcf, 0);
        }
        else
        {
            amrex::Vector<amrex::Box> cloc_boxes;
            for (amrex::MFIter mfi(mf[p->level-1]); mfi.isValid(); ++mfi)
            {
                cloc_boxes.push_back(mfi.validbox());
            }

            amrex::Gpu::DeviceVector<amrex::Box> cdevice_boxes(cloc_boxes.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, cloc_boxes.begin(), cloc_boxes.end(), cdevice_boxes.begin());

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> cbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{face_bc_values, face_heat_values, gcv, p->margin, static_cast<bool>(p->j_dir), bc_decision,
                                                        cdevice_boxes.data(), static_cast<int>(cdevice_boxes.size())});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> cphysbcf(
                p->amrex_geometry[p->level-1], BCRecs[p->level-1], cbf);

            amrex::Vector<amrex::Box> floc_boxes;
            for (amrex::MFIter mfi(mf[p->level]); mfi.isValid(); ++mfi)
            {
                floc_boxes.push_back(mfi.validbox());
            }

            amrex::Gpu::DeviceVector<amrex::Box> fdevice_boxes(floc_boxes.size());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, floc_boxes.begin(), floc_boxes.end(), fdevice_boxes.begin());

            amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>> fbf(
                amrex_bc_func::MyExtBCFillField<BCDecision>{face_bc_values, face_heat_values, gcv, p->margin, static_cast<bool>(p->j_dir), bc_decision,
                                                        fdevice_boxes.data(), static_cast<int>(fdevice_boxes.size())});

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<amrex_bc_func::MyExtBCFillField<BCDecision>>> fphysbcf(
                p->amrex_geometry[p->level], BCRecs[p->level], fbf);

            amrex::CellConservativeLinear mapper;
            const amrex::IntVect ratio = 2 * amrex::IntVect::TheUnitVector();

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
}

#endif
