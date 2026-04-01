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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef LOOPING_H_
#define LOOPING_H_

#include"definitions.h"

#include"looping2D.h"

#include"iterators1D.h"
#include"iterators2D.h"
#include"iterators3D.h"
#if USE_AMREX
    #include <AMReX_MFIter.H>
    #include <AMReX_MultiFab.H>
    #include <AMReX_GpuLaunchFunctsC.H>
#endif

// LOOPs
#if USE_AMREX
    // -------------------------------------------------------------------------
    // FIELDLOOP — iterate all AMReX levels over a mutable field, expose const
    // fields as named MultiFab refs and Array4 locals (both at tile scope),
    // then launch a ParallelFor with the given body.
    //
    //   mut          – mutable field; shadowed by an Array4 of the same name so
    //                  'mut(i,j,k)' writes inside the lambda.
    //   const_decls  – semicolon-separated FIELD_CONST / FIELD_CONST_MEMBER calls.
    //   body         – expression evaluated per (i,j,k); must end with ';'.
    //
    // Helpers for const_decls:
    //   FIELD_CONST(name)            – simple field; creates _fl_mf_name ref and
    //                                  shadows 'name' with a const Array4.
    //   FIELD_CONST_MEMBER(ptr,mbr)  – ptr->mbr field; creates _fl_mf_a_mbr ref
    //                                  and 'member_mbr' const Array4; body uses member_mbr(i,j,k).
    //
    // Example:
    //   FIELDLOOP(ark2,
    //       FIELD_CONST(ls); FIELD_CONST(ark1); FIELD_CONST_MEMBER(a, L),
    //       ark2(i,j,k) = 0.75*ls(i,j,k) + 0.25*ark1(i,j,k) + 0.25*dt*a_L(i,j,k);)
    // -------------------------------------------------------------------------
    #ifdef AMREX_USE_OMP
    #  define _FL_OMP _Pragma("omp parallel if (amrex::Gpu::notInLaunchRegion())")
    #else
    #  define _FL_OMP
    #endif

    #define FIELD_CONST(name) \
        const auto& _fl_mf_##name = (name).GetMultiFab(_fl_lev); \
        const auto& name = _fl_mf_##name[_fl_mfi].const_array();

    #define FIELD_CONST_MEMBER(ptr, member) \
        const auto& _fl_mf_a_##member = (ptr)->member.GetMultiFab(_fl_lev); \
        const auto& member_##member = _fl_mf_a_##member[_fl_mfi].const_array();

    #define _FIELDLOOP_IMPL(mut_expr, mut_name, const_decls, body) \
        for (int _fl_lev = p->nlevs - 1; _fl_lev >= 0; --_fl_lev) \
        { \
            auto& _fl_mf = (mut_expr).GetMultiFab(_fl_lev); \
            _FL_OMP \
            for (amrex::MFIter _fl_mfi(_fl_mf, amrex::TilingIfNotGPU()); \
                 _fl_mfi.isValid(); ++_fl_mfi) \
            { \
                const amrex::Box& _fl_bx = _fl_mfi.tilebox(); \
                auto const& mut_name = _fl_mf[_fl_mfi].array(); \
                const_decls; \
                amrex::ParallelFor(_fl_bx, \
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) { body }); \
            } \
        }

    #define FIELDLOOP(mut, const_decls, body) \
        _FIELDLOOP_IMPL(mut, mut, const_decls, body)

    #define FIELDLOOP_MEMBER(ptr, member, const_decls, body) \
        _FIELDLOOP_IMPL((ptr)->member, member, const_decls, body)

    // -------------------------------------------------------------------------
    // FIELDLOOP_INC — like FIELDLOOP but also maintains increment::i/j/k inside
    // the lambda and exposes const fields as LocalArr4Const (offset-aware) so
    // that legacy helpers using increment::i/j/k work on GPU/CPU uniformly.
    //
    // Warning: mut is expecting access using global indices (i,j,k)
    //   since it is not wrapped in a LocalArr4.
    //
    // Differences from FIELDLOOP:
    //   • sets p->level = lev per level (and resets to 0 after the loop)
    //   • calls p->set_tile_mfi per tile (and restores default after the loop)
    //   • computes ox/oy/oz = tilebox.smallEnd for offset subtraction
    //   • const_decls should use FIELD_CONST_INC / FIELD_CONST_MEMBER_INC
    //   • lambda sets increment::i/j/k = i-ox / j-oy / k-oz before body
    //
    // Helpers for const_decls (require field.h to be included at call site):
    //   FIELD_CONST_INC(name)            – simple field → LocalArr4Const 'name'
    //   FIELD_CONST_MEMBER_INC(ptr,mbr)  – ptr->mbr    → LocalArr4Const 'a_mbr'
    //
    // Example:
    //   FIELDLOOP_INC(a->L,
    //       FIELD_CONST_INC(b); FIELD_CONST_INC(uvel); FIELD_CONST_MEMBER_INC(a, H),
    //       L(i,j,k) += aij(p, a, b, 4, uvel, a_H, ...);)
    // -------------------------------------------------------------------------
    #define FIELD_CONST_INC(name) \
        const auto& _fl_mf_##name = (name).GetMultiFab(_fl_lev); \
        const LocalArr4Const name(_fl_mf_##name[_fl_mfi].const_array(), ox, oy, oz);

    #define FIELD_CONST_MEMBER_INC(ptr, member) \
        const auto& _fl_mf_a_##member = (ptr)->member.GetMultiFab(_fl_lev); \
        const LocalArr4Const member_##member(_fl_mf_a_##member[_fl_mfi].const_array(), ox, oy, oz);

    #define _FIELDLOOP_INC_IMPL(mut_expr, mut_name, const_decls, body) \
        for (int _fl_lev = p->nlevs - 1; _fl_lev >= 0; --_fl_lev) \
        { \
            p->level = _fl_lev; \
            auto& _fl_mf = (mut_expr).GetMultiFab(_fl_lev); \
            _FL_OMP \
            for (amrex::MFIter _fl_mfi(_fl_mf, amrex::TilingIfNotGPU()); \
                 _fl_mfi.isValid(); ++_fl_mfi) \
            { \
                const amrex::Box& _fl_bx = _fl_mfi.tilebox(); \
                const int ox = _fl_bx.smallEnd(0); \
                const int oy = _fl_bx.smallEnd(1); \
                const int oz = _fl_bx.smallEnd(2); \
                p->set_tile_mfi(&_fl_mfi); \
                auto const& mut_name = _fl_mf[_fl_mfi].array(); \
                const_decls; \
                amrex::ParallelFor(_fl_bx, \
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) { \
                        increment::i = i - ox; \
                        increment::j = j - oy; \
                        increment::k = k - oz; \
                        body \
                    }); \
            } \
        } \
        p->level = 0; \
        p->set_tile_mfi(p->default_cell_mfi.get());

    #define FIELDLOOP_INC(mut, const_decls, body) \
        _FIELDLOOP_INC_IMPL(mut, mut, const_decls, body)

    #define FIELDLOOP_INC_MEMBER(ptr, member, const_decls, body) \
        _FIELDLOOP_INC_IMPL((ptr)->member, member, const_decls, body)

    #define LEVEL_LOOP \
        for (struct { lexer* ctx; bool active; } \
                _level_guard{p, true}; \
            _level_guard.active; \
            _level_guard.ctx->level = 0, _level_guard.active = false) \
            for (_level_guard.ctx->level = 0; \
                _level_guard.ctx->level < _level_guard.ctx->nlevs; \
                ++_level_guard.ctx->level)
    #define TILE_LOOP \
        for (amrex::MFIter _tile_mfi(p->amr_cell_mf[p->level],amrex::TilingIfNotGPU()); _tile_mfi.isValid(); ++_tile_mfi) \
            for (struct { lexer* ctx; amrex::MFIter* saved; } \
                    _guard{p, p->set_tile_mfi(&_tile_mfi)}; \
                _guard.ctx != nullptr; \
                _guard.ctx->set_tile_mfi(_guard.saved ? _guard.saved : _guard.ctx->default_cell_mfi.get()), \
                _guard.ctx = nullptr)

    #define IMAX_LOOP (p->amr_tile_hi.x - p->amr_tile_lo.x)
    #define JMAX_LOOP (p->amr_tile_hi.y - p->amr_tile_lo.y)
    #define KMAX_LOOP (p->amr_tile_hi.z - p->amr_tile_lo.z)
    #define MARGIN_I p->amr_cell_mf[p->level].nGrow(0)
    #define MARGIN_J p->amr_cell_mf[p->level].nGrow(1)
    #define MARGIN_K p->amr_cell_mf[p->level].nGrow(2)
#else
    // Non-AMReX FIELDLOOP: simple fields are used directly from the body.
    // ptr->member fields get a local alias so a_member(i,j,k) works in both modes.
    #define FIELD_CONST(name)
    #define FIELD_CONST_MEMBER(ptr, member)           auto& member_##member = (ptr)->member;
    #define FIELDLOOP(mut, const_decls, body)         { const_decls; PLAINLOOP PCHECK body }
    #define FIELDLOOP_MEMBER(ptr, member, const_decls, body) \
        { auto& member = (ptr)->member; const_decls; PLAINLOOP PCHECK body }

    #define FIELD_CONST_INC(name)
    #define FIELD_CONST_MEMBER_INC(ptr, member)       auto& member_##member = (ptr)->member;
    #define FIELDLOOP_INC(mut, const_decls, body)     { const_decls; PLAINLOOP PCHECK body }
    #define FIELDLOOP_INC_MEMBER(ptr, member, const_decls, body) \
        { auto& member = (ptr)->member; const_decls; PLAINLOOP PCHECK body }

    #define LEVEL_LOOP
    #define TILE_LOOP
    #define IMAX_LOOP p->knox-1
    #define JMAX_LOOP p->knoy-1
    #define KMAX_LOOP p->knoz-1
    #define MARGIN_I p->margin
    #define MARGIN_J p->margin
    #define MARGIN_K p->margin
#endif

#define ILOOP for (i = 0; i <= IMAX_LOOP; ++i)
#define JLOOP for (j = 0; j <= JMAX_LOOP; ++j)
#define KLOOP for (k = 0; k <= KMAX_LOOP; ++k)

#define ITLOOP for(i = 0; i <= IMAX_LOOP+1; ++i)
#define JTLOOP for(j = 0; j <= JMAX_LOOP+1; ++j)
#define KTLOOP for(k = 0; k <= KMAX_LOOP+1; ++k)

#define ITPLOOP for(i = -1; i <= IMAX_LOOP; ++i)
#define JTPLOOP for(j = -1; j <= JMAX_LOOP; ++j)
#define KTPLOOP for(k = -1; k <= KMAX_LOOP; ++k)

#define IBLOOP for(i = -1; i <= IMAX_LOOP+1; ++i)
#define JBLOOP for(j = -1; j <= JMAX_LOOP+1; ++j)
#define KBLOOP for(k = -1; k <= KMAX_LOOP+1; ++k)

#define IMALOOP for(i = -MARGIN_I; i <= IMAX_LOOP + MARGIN_I; ++i)
#define JMALOOP for(j = -MARGIN_J; j <= JMAX_LOOP + MARGIN_J; ++j)
#define KMALOOP for(k = -MARGIN_K; k <= KMAX_LOOP + MARGIN_K; ++k)

#define IFLEXLOOP for(i = 0; i <= IMAX_LOOP - ulast; ++i)
#define JFLEXLOOP for(j = 0; j <= JMAX_LOOP - vlast; ++j)
#define KFLEXLOOP for(k = 0; k <= KMAX_LOOP - wlast; ++k)

#define IULOOP for(i = 0; i <= IMAX_LOOP - p->ulast; ++i)
#define JVLOOP for(j = 0; j <= JMAX_LOOP - p->vlast; ++j)
#define KWLOOP for(k = 0; k <= KMAX_LOOP - p->wlast; ++k)

#define FILOOP ILOOP
#define FJLOOP JLOOP
#define FKLOOP for(k = 0; k <= KMAX_LOOP + 1; ++k)

#define ETALOC  for(k=a->etaloc(i,j); k<a->etaloc(i,j)+1; ++k)
#define FETALOC for(k=c->etaloc(i,j); k<c->etaloc(i,j)+1; ++k)

// CONDITIONS
#define FLEXCHECK   if(flag[IJK]>0)
#if USE_AMREX
    // #define UCHECK      if(p->flag1(i,j,k)>0)
    // #define UFLUIDCHECK if(p->flag1(i,j,k)>=AIR_FLAG)
    // #define USCHECK     if(p->flag1(i,j,k)<0)
    // #define VCHECK      if(p->flag2(i,j,k)>0)
    // #define VFLUIDCHECK if(p->flag2(i,j,k)>=AIR_FLAG)
    // #define VSCHECK     if(p->flag2(i,j,k)<0)
    // #define WCHECK      if(p->flag3(i,j,k)>0)
    // #define WFLUIDCHECK if(p->flag3(i,j,k)>=AIR_FLAG)
    // #define WSCHECK     if(p->flag3(i,j,k)<0)
    // #define PCHECK      if(p->flag4(i,j,k)>0)
    // #define SCHECK      if(p->flag4(i,j,k)<0)
    // #define PFLUIDCHECK if(p->flag4(i,j,k)>=AIR_FLAG)
    // #define PAIR_CHECK  if(p->flag4(i,j,k)==AIR_FLAG)
    // #define SFLUIDCHECK if(p->flag4(i,j,k)<AIR_FLAG)
    // #define PSOLIDCHECK if(p->flag4(i,j,k)>SOLID_FLAG)
    // #define PBASECHECK  if(p->flag4(i,j,k)>OBJ_FLAG)
    #define UCHECK
    #define UFLUIDCHECK
    #define USCHECK
    #define VCHECK
    #define VFLUIDCHECK
    #define VSCHECK
    #define WCHECK
    #define WFLUIDCHECK
    #define WSCHECK
    #define PCHECK
    #define SCHECK
    #define PFLUIDCHECK
    #define PAIR_CHECK
    #define SFLUIDCHECK
    #define PSOLIDCHECK
    #define PBASECHECK
#else
    #define UCHECK      if(p->flag1[IJK]>0)
    #define UFLUIDCHECK if(p->flag1[IJK]>=AIR_FLAG)
    #define USCHECK     if(p->flag1[IJK]<0)
    #define VCHECK      if(p->flag2[IJK]>0)
    #define VFLUIDCHECK if(p->flag2[IJK]>=AIR_FLAG)
    #define VSCHECK     if(p->flag2[IJK]<0)
    #define WCHECK      if(p->flag3[IJK]>0)
    #define WFLUIDCHECK if(p->flag3[IJK]>=AIR_FLAG)
    #define WSCHECK     if(p->flag3[IJK]<0)
    #define PCHECK      if(p->flag4[IJK]>0)
    #define SCHECK      if(p->flag4[IJK]<0)
    #define PFLUIDCHECK if(p->flag4[IJK]>=AIR_FLAG)
    #define PAIR_CHECK  if(p->flag4[IJK]==AIR_FLAG)
    #define SFLUIDCHECK if(p->flag4[IJK]<AIR_FLAG)
    #define PSOLIDCHECK if(p->flag4[IJK]>SOLID_FLAG)
    #define PBASECHECK  if(p->flag4[IJK]>OBJ_FLAG)
#endif
#define FPCHECK     if(p->flag7[FIJK]>0)
#define FPWDCHECK   if(p->flag7[FIJK]>0 && p->wet[IJ]>0)
#define FSCHECK     if(p->flag7[FIJK]<=0)
#define FSWDCHECK   if(p->flag7[FIJK]<=0 || p->wet[IJ]==0)

// POROSITY
#define PORVAL1 (0.5*(a->porosity(i+1,j,k) + a->porosity(i,j,k)))
#define PORVAL2 (0.5*(a->porosity(i,j+1,k) + a->porosity(i,j,k)))
#define PORVAL3 (0.5*(a->porosity(i,j,k+1) + a->porosity(i,j,k)))
#define PORVAL4 a->porosity(i,j,k)

#define PORVAL4px a->porosity(i+1,j,k)
#define PORVAL4py a->porosity(i,j+1,k)
#define PORVAL4pz a->porosity(i,j,k+1)

#define CPOR4px   (1.0/(1.0+(p->B260*(PORVAL4px<1.0?1.0:0.0))))
#define CPOR4py   (1.0/(1.0+(p->B260*(PORVAL4py<1.0?1.0:0.0))))
#define CPOR4pz   (1.0/(1.0+(p->B260*(PORVAL4pz<1.0?1.0:0.0))))

#define CPOR1   (1.0/(1.0+(p->B260*(PORVAL1<1.0?1.0:0.0))))
#define CPOR2   (1.0/(1.0+(p->B260*(PORVAL2<1.0?1.0:0.0))))
#define CPOR3   (1.0/(1.0+(p->B260*(PORVAL3<1.0?1.0:0.0))))
#define CPOR4   (1.0/(1.0+(p->B260*(PORVAL4<1.0?1.0:0.0))))

#define PORVAL1m (0.5*(a->porosity(i,j,k) + a->porosity(i-1,j,k)))
#define PORVAL2m (0.5*(a->porosity(i,j,k) + a->porosity(i,j-1,k)))
#define PORVAL3m (0.5*(a->porosity(i,j,k) + a->porosity(i,j,k-1)))

#define CPOR1m   (1.0/(1.0+(p->B260*(PORVAL1m<1.0?1.0:0.0))))
#define CPOR2m   (1.0/(1.0+(p->B260*(PORVAL2m<1.0?1.0:0.0))))
#define CPOR3m   (1.0/(1.0+(p->B260*(PORVAL3m<1.0?1.0:0.0))))

#define PORVAL1p (0.5*(a->porosity(i+2,j,k) + a->porosity(i+1,j,k)))
#define PORVAL2p (0.5*(a->porosity(i,j+2,k) + a->porosity(i,j+1,k)))
#define PORVAL3p (0.5*(a->porosity(i,j,k+2) + a->porosity(i,j,k+1)))

#define CPOR1p   (1.0/(1.0+(p->B260*(PORVAL1p<1.0?1.0:0.0))))
#define CPOR2p   (1.0/(1.0+(p->B260*(PORVAL2p<1.0?1.0:0.0))))
#define CPOR3p   (1.0/(1.0+(p->B260*(PORVAL3p<1.0?1.0:0.0))))

// COMBINDED LOOPS
#define MultiGridLOOP LEVEL_LOOP TILE_LOOP
#define IJKLOOP ILOOP JLOOP KLOOP
#define KJILOOP KLOOP JLOOP ILOOP

// MAIN LOOPS
#define PLAINLOOP MultiGridLOOP IJKLOOP
#define LOOP PLAINLOOP PCHECK
#define BASELOOP PLAINLOOP PBASECHECK
#define BASEREVLOOP MultiGridLOOP KJILOOP PBASECHECK
#define TPLOOP KTPLOOP JTPLOOP ITPLOOP

#define KYLREVLOOP for(k=p->knoz-1; k>=1; --k)
#define YLLOOP ILOOP JLOOP KYLREVLOOP PBASECHECK

#define ALOOP PLAINLOOP PSOLIDCHECK
#define AIRLOOP PLAINLOOP PAIR_CHECK

#define MALOOP IMALOOP JMALOOP KMALOOP

// BOUNDARY LOOPS
#define BLOOP MultiGridLOOP IBLOOP JBLOOP KBLOOP
#define BBASELOOP MultiGridLOOP IBLOOP JBLOOP KBLOOP PBASECHECK

// FLUID LOOPS
#define ULOOP MultiGridLOOP IULOOP JLOOP KLOOP UCHECK
#define VLOOP MultiGridLOOP ILOOP JVLOOP KLOOP VCHECK
#define WLOOP MultiGridLOOP ILOOP JLOOP KWLOOP WCHECK

#define UFLUIDLOOP MultiGridLOOP IULOOP JLOOP KLOOP UFLUIDCHECK
#define VFLUIDLOOP MultiGridLOOP ILOOP JVLOOP KLOOP VFLUIDCHECK
#define WFLUIDLOOP MultiGridLOOP ILOOP JLOOP KWLOOP WFLUIDCHECK
#define FLUIDLOOP PLAINLOOP PFLUIDCHECK

// SOLVER LOOPS
#define FLEXLOOP MultiGridLOOP IFLEXLOOP JFLEXLOOP KFLEXLOOP FLEXCHECK

// FNPF LOOPS
#define FKJILOOP MultiGridLOOP FKLOOP JLOOP ILOOP
#define FLOOP MultiGridLOOP ILOOP JLOOP FKLOOP FPCHECK
#define FBASELOOP MultiGridLOOP ILOOP JLOOP FKLOOP
#define FILOOP4 MultiGridLOOP ILOOP JLOOP ETALOC PFLUIDCHECK
#define FFILOOP4 MultiGridLOOP for(k=KMAX_LOOP+1; k<KMAX_LOOP+2; ++k) ILOOP JLOOP PSLICECHECK4

// FORCE LOOP
#define NDBASELOOP MultiGridLOOP ITPLOOP JTPLOOP KTPLOOP

#define NETLOOP for (int n=0; n<p->net_count; ++n)

#define NLOOP1 for(n=p->sizeM1[0]; n<p->sizeM1[1]; ++n)
#define NLOOP2 for(n=p->sizeM2[0]; n<p->sizeM2[1]; ++n)
#define NLOOP3 for(n=p->sizeM3[0]; n<p->sizeM3[1]; ++n)
#define NLOOP4 for(n=p->sizeM4[0]; n<p->sizeM4[1]; ++n)
#define NLOOP4A for(n=p->sizeM4a[0]; n<p->sizeM4a[1]; ++n)
#define NLOOP6 for(n=p->sizeM6[0]; n<p->sizeM6[1]; ++n)
#define NLOOP9 for(n=p->sizeM9[0]; n<p->sizeM9[1]; ++n)
#define NLOOP for(n=sizeM[0]; n<sizeM[1]; ++n)

//MAX, MIN, SIGN
#define MAX(aAa,bBb) ((aAa)>(bBb)?(aAa):(bBb))
#define MIN(aAa,bBb) ((aAa)<(bBb)?(aAa):(bBb))
#define SIGN(aAa) ((aAa)>= 0.0 ? 1.0 : -1.0)

#define MAX3(aAa,bBb,cCc) (((aAa)>(bBb)?(aAa):(bBb))>cCc?((aAa)>(bBb)?(aAa):(bBb)):cCc)
#define MIN3(aAa,bBb,cCc) (((aAa)<(bBb)?(aAa):(bBb))<cCc?((aAa)<(bBb)?(aAa):(bBb)):cCc)

//GCB
#define GCB1 for(n=0;n<p->gcb1_count;++n)
#define GCB1CHECK if(p->gcb1[n][3]>0)
#define GC1LOOP GCB1 GCB1CHECK

#define QGCB1 for(q=0;q<p->gcb1_count;++q)
#define QGCB1CHECK if(p->gcb1[q][3]>0)
#define QGC1LOOP QGCB1 QGCB1CHECK

#define QQGCB1 for(qq=0;qq<p->gcb1_count;++qq)
#define QQGCB1CHECK if(p->gcb1[qq][3]>0)
#define QQGC1LOOP QQGCB1 QQGCB1CHECK

#define GCB2 for(n=0;n<p->gcb2_count;++n)
#define GCB2CHECK if(p->gcb2[n][3]>0)
#define GC2LOOP GCB2 GCB2CHECK

#define QGCB2 for(q=0;q<p->gcb2_count;++q)
#define QGCB2CHECK if(p->gcb2[q][3]>0)
#define QGC2LOOP QGCB2 QGCB2CHECK

#define QQGCB2 for(qq=0;qq<p->gcb2_count;++qq)
#define QQGCB2CHECK if(p->gcb2[qq][3]>0)
#define QQGC2LOOP QQGCB2 QQGCB2CHECK

#define GCB3 for(n=0;n<p->gcb3_count;++n)
#define GCB3CHECK if(p->gcb3[n][3]>0)
#define GC3LOOP GCB3 GCB3CHECK

#define QGCB3 for(q=0;q<p->gcb3_count;++q)
#define QGCB3CHECK if(p->gcb3[q][3]>0)
#define QGC3LOOP QGCB3 QGCB3CHECK

#define QQGCB3 for(qq=0;qq<p->gcb3_count;++qq)
#define QQGCB3CHECK if(p->gcb3[qq][3]>0)
#define QQGC3LOOP QQGCB3 QQGCB3CHECK

#define GCB4 for(n=0;n<p->gcb4_count;++n)
#define GCB4CHECK if(p->gcb4[n][3]>0)
#define GC4LOOP GCB4 GCB4CHECK

#define QGCB4 for(q=0;q<p->gcb4_count;++q)
#define QGCB4CHECK if(p->gcb4[q][3]>0)
#define QGC4LOOP QGCB4 QGCB4CHECK

#define QQGCB4 for(qq=0;qq<p->gcb4_count;++qq)
#define QQGCB4CHECK if(p->gcb4[qq][3]>0)
#define QQGC4LOOP QQGCB4 QQGCB4CHECK

//df
#define QGCDF1 for(q=0;q<p->gcdf1_count;++q)
#define QGCDF1CHECK if(p->gcdf1[q][3]>0)
#define QGCDF1LOOP QGCDF1 QGCDF1CHECK

#define GCDF1 for(n=0;n<p->gcdf1_count;++n)
#define GCDF1CHECK if(p->gcdf1[n][3]>0)
#define GCDF1LOOP GCDF1 GCDF1CHECK

#define QGCDF2 for(q=0;q<p->gcdf2_count;++q)
#define QGCDF2CHECK if(p->gcdf2[q][3]>0)
#define QGCDF2LOOP QGCDF2 QGCDF2CHECK

#define GCDF2 for(n=0;n<p->gcdf2_count;++n)
#define GCDF2CHECK if(p->gcdf2[n][3]>0)
#define GCDF2LOOP GCDF2 GCDF2CHECK

#define QGCDF3 for(q=0;q<p->gcdf3_count;++q)
#define QGCDF3CHECK if(p->gcdf3[q][3]>0)
#define QGCDF3LOOP QGCDF3 QGCDF3CHECK

#define GCDF3 for(n=0;n<p->gcdf3_count;++n)
#define GCDF3CHECK if(p->gcdf3[n][3]>0)
#define GCDF3LOOP GCDF3 GCDF3CHECK

#define QGCDF4 for(q=0;q<p->gcdf4_count;++q)
#define QGCDF4CHECK if(p->gcdf4[q][3]>0)
#define QGCDF4LOOP QGCDF4 QGCDF4CHECK

#define GCDF4 for(n=0;n<p->gcdf4_count;++n)
#define GCDF4CHECK if(p->gcdf4[n][3]>0)
#define GCDF4LOOP GCDF4 GCDF4CHECK

#endif
