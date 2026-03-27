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

#include"iterators.h"
#include"definitions.h"
#include"looping2D.h"
#include"iterators2D.h"
#include"iterators1D.h"
#if USE_AMREX
    #include <AMReX_MFIter.H>
    #include <AMReX_MultiFab.H>

    // LOOPs
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
#define FPCHECK     if(p->flag7[FIJK]>0)
#define FPWDCHECK   if(p->flag7[FIJK]>0 && p->wet[IJ]>0)
#define FSCHECK     if(p->flag7[FIJK]<=0)
#define FSWDCHECK   if(p->flag7[FIJK]<=0 || p->wet[IJ]==0)

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
#define FFILOOP4 MultiGridLOOP ILOOP JLOOP FETALOC FPCHECK

// FORCE LOOP
#define NDBASELOOP MultiGridLOOP ITPLOOP JTPLOOP KTPLOOP

#define NETLOOP for (int n=0; n<p->net_count; ++n)

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
