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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#if USE_AMREX
#ifndef FIELDS_AMReX_H_
#define FIELDS_AMReX_H_

#include "field_amrex.h"

class field1 : public field_amrex
{
public:
    field1(lexer*);
    field1(lexer*, amrex::Vector<amrex::MultiFab>* shared_mf, int comp);
    virtual ~field1() = default;
    void FillDomainBoundary(int gcv) override;
    void UpdateBCRecs(int gcv) override;
private:
    amrex_bc_func::Field1BcDecision::Field1Params params;
    void init_params(lexer*);
};

class field2 : public field_amrex
{
public:
    field2(lexer*);
    field2(lexer*, amrex::Vector<amrex::MultiFab>* shared_mf, int comp);
    virtual ~field2() = default;
    void FillDomainBoundary(int gcv) override;
    void UpdateBCRecs(int gcv) override;
private:
    amrex_bc_func::Field2BcDecision::Field2Params params;
    void init_params(lexer*);
};

class field3 : public field_amrex
{
public:
    field3(lexer*);
    field3(lexer*, amrex::Vector<amrex::MultiFab>* shared_mf, int comp);
    virtual ~field3() = default;
    void FillDomainBoundary(int gcv) override;
    void UpdateBCRecs(int gcv) override;
private:
    amrex_bc_func::Field3BcDecision::Field3Params params;
    void init_params(lexer*);
};

class field4 : public field_amrex
{
public:
    field4(lexer*);
    field4(lexer*, amrex::Vector<amrex::MultiFab>* shared_mf, int comp);
    virtual ~field4() = default;
    void FillDomainBoundary(int gcv) override;
    void UpdateBCRecs(int gcv) override;
private:
    amrex_bc_func::Field4BcDecision::Field4Params params;
    void init_params(lexer*);
};

/*!
 * @brief Sigma-grid vertical-node field: n*m*(o+1) planes.
 *
 * AMReX counterpart of the flat-array field7. Backed by a z-nodal MultiFab --
 * amrex::convert(amrex_box_array[lev], IntVect(0,0,1)), selected by
 * DataLocation::NODE_Z -- which shares box count, ordering and
 * DistributionMap with every cell-centred field. That is what lets it be
 * driven by the same MFIter and the same TileCtx: field access resolves
 * through get_array(level, MFIter::LocalIndex()), and LocalIndex is a pure
 * function of the DistributionMap.
 *
 * The vertical index needs no translation either. TILE_LOOP sets
 * m_cached_oz = amr_tile_lo.z, and amrex::convert preserves smallEnd, so
 * FKLOOP's local k = 0..KMAX_LOOP+1 maps onto the box's valid nodal range
 * loz..loz+nz exactly -- one plane more than KLOOP walks, which is the point.
 *
 * The cost of that overlap: the plane on a z-split box seam is valid in TWO
 * boxes. field_amrex::FillBoundary OverrideSyncs NODE_Z fields before filling
 * ghosts so the duplicates cannot disagree.
 */
class field7 : public field_amrex
{
public:
    field7(lexer*);
    field7(lexer*, amrex::Vector<amrex::MultiFab>* shared_mf, int comp);
    virtual ~field7() = default;
    void FillDomainBoundary(int gcv) override;
    void UpdateBCRecs(int gcv) override;
private:
    /// Field7BcDecision aliases Field4BcDecision: a sigma-grid scalar's
    /// domain-boundary rules are field4's. See amrex_bc_func.h.
    amrex_bc_func::Field7BcDecision::Field4Params params;
    void init_params(lexer*);
};

#endif
#endif