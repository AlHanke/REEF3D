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
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>

class amrex_bc_func
{
public:
    amrex_bc_func() = default;
    virtual ~amrex_bc_func() = default;

    struct MyExtBCFill {
        AMREX_GPU_HOST_DEVICE
        MyExtBCFill() = default;

        AMREX_GPU_HOST_DEVICE
        MyExtBCFill(const amrex::Array<int,6>& values, const int gcv, const amrex::Box* boxes, int num_boxes)
            : bc_values(values), m_gcv(gcv), m_boxes(boxes), m_num_boxes(num_boxes) {}

        AMREX_GPU_DEVICE
        void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                        const int dcomp, const int numcomp,
                        amrex::GeometryData const& geom, const amrex::Real time,
                        const amrex::BCRec* bcr, const int bcomp,
                        const int orig_comp) const
        {
            amrex::ignore_unused(time, orig_comp);

            bool should_fill = false;
            int face_for_bc = 0;

            // Check against local boxes
            for (int i = 0; i < m_num_boxes; ++i)
            {
                const amrex::Box& box = m_boxes[i];
                int f = detect_face(iv, box);
                if (f > 0)
                {
                    // Face projection of a local box -> Fill
                    should_fill = true;
                    face_for_bc = f;
                    break;
                }
                else if (is_corner_layer1(iv, box))
                {
                    // Corner of a local box, layer 1 -> Fill
                    should_fill = true;
                }
            }

            amrex::IntVect interior = iv;
            const amrex::Box dom = geom.Domain();
            for(int dir=0; dir<AMREX_SPACEDIM; ++dir)
            {
                if(interior[dir] < dom.smallEnd(dir))
                    interior[dir] = dom.smallEnd(dir);
                else if(interior[dir] > dom.bigEnd(dir))
                    interior[dir] = dom.bigEnd(dir);
            }

            if (should_fill)
            {
                bool bc_active = true;
                if (face_for_bc > 0) {
                    if (bc_values[face_for_bc-1] == 0) bc_active = false;
                }

                if (bc_active)
                {
                    for(int n=0; n<numcomp; ++n)
                    {
                        const amrex::BCRec& bc = bcr[bcomp + n];
                        amrex::ignore_unused(bc);

                        // Neumann
                        dest(iv, dcomp+n) = dest(interior, dcomp+n);
                    }
                }
            }
        }

    private:
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

        amrex::Array<int,2*AMREX_SPACEDIM> bc_values{};
        const amrex::Box* m_boxes = nullptr;
        int m_num_boxes = 0;
        int m_gcv = 0;
    };
};

#endif
#endif
