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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef LOOPING2D_H_
#define LOOPING2D_H_

// SLICE CONDITIONS
#if USE_AMREX
    #define PSLICECHECK1  if(p->flagslice1(i,j)>0)
    #define PSLICECHECK2  if(p->flagslice2(i,j)>0)
    #define PSLICECHECK4  if(p->flagslice4(i,j)>0)
    #define SSLICECHECK4  if(p->flagslice4(i,j)<0)
#else
    #define PSLICECHECK1  if(p->flagslice1[IJ]>0)
    #define PSLICECHECK2  if(p->flagslice2[IJ]>0)
    #define PSLICECHECK4  if(p->flagslice4[IJ]>0)
    #define SSLICECHECK4  if(p->flagslice4[IJ]<0)
#endif

#define WETDRYCHK if(p->wet[IJ]>0)
#define SEDSLICECHECK if(s->DFBED[IJ]>0)
#define SLICEFLEXCHECK  if(flagslice[IJ]>0)

#define WETDRY1 if(b->wet1(i,j)==1)
#define WETDRY2 if(b->wet2(i,j)==1)
#define WETDRYDEEP1 if(b->wet1(i,j)==1 && b->deep1(i,j)==1)
#define WETDRYDEEP2 if(b->wet2(i,j)==1 && b->deep2(i,j)==1)
#define WETDRY if(p->wet[IJ]==1)
#define WETDRYDEEP if(p->wet[IJ]==1 && p->deep[IJ]==1)
#define WETDRYDEEPBREAK if(p->wet[IJ]==1 && p->deep[IJ]==1 && d->breaking(i,j)==0)

// SLICE BASE LOOPS
#define SLICEBASELOOP MultiGridLOOP ILOOP JLOOP
#define JILOOP JLOOP ILOOP

// SLICE LOOPS
#define SLICELOOP1 MultiGridLOOP IULOOP JLOOP  PSLICECHECK1
#define SLICELOOP2 MultiGridLOOP ILOOP JVLOOP  PSLICECHECK2
#define SLICELOOP4 SLICEBASELOOP  PSLICECHECK4

#define SEDSLICELOOP SLICEBASELOOP  PSLICECHECK4 SEDSLICECHECK

// SLICE OWNER LOOPS — visit each (i,j) column of the level EXACTLY ONCE.
//
// The slice loops above do not: slice data lives on z-slabs of the 3D BoxArray,
// so every box in a column flattens onto the same (i,j), and MultiGridLOOP walks
// them all. With MFIter_TILING live it is worse than one visit per box — the
// default mfiter_tile_size splits a box every 8 cells in z, and since these loops
// never iterate k, each z-tile repeats the whole (i,j) footprint of its box.
//
// For a loop that WRITES a column value that is the right behaviour: the view is
// overlapping and is only re-broadcast at the next makeUnique, so every copy has
// to be written or a reader on a non-owner tile sees stale data. Use SLICELOOP*
// there — this is not a stricter drop-in replacement.
//
// Use SLICEOWNLOOP* wherever a repeat is wrong rather than redundant: emitting one
// record per column (mgcslice*::gcb_seed), counting, or accumulating with +=.
//
// PSLICEOWNER sits between the tile loop and the (i,j) sweep so a non-owner tile
// is skipped whole. See TileCtx::slice_owner for the predicate, and
// grid_amrex::build_slice_owner for the box-granular mask it agrees with.
#if USE_AMREX
    #define PSLICEOWNER if(p->amr_slice_owner)
#else
    #define PSLICEOWNER
#endif

#define SLICEOWNBASELOOP MultiGridLOOP PSLICEOWNER ILOOP JLOOP

#define SLICEOWNLOOP1 MultiGridLOOP PSLICEOWNER IULOOP JLOOP  PSLICECHECK1
#define SLICEOWNLOOP2 MultiGridLOOP PSLICEOWNER ILOOP  JVLOOP PSLICECHECK2
#define SLICEOWNLOOP4 SLICEOWNBASELOOP  PSLICECHECK4

#define TPSLICELOOP ITPLOOP JTPLOOP

#define SLICEFLEXLOOP IFLEXLOOP JFLEXLOOP SLICEFLEXCHECK

// GCBSL

#define GCSLB1 for(n=0;n<p->gcbsl1.ssize(p->level);++n)
#define GCSLB1CHECK if(p->gcbsl1[p->level][n].cs>0)
#define GCSL1LOOP GCSLB1 GCSLB1CHECK

#define QGCSLB1 for(q=0;q<p->gcbsl1.ssize(p->level);++q)
#define QGCSLB1CHECK if(p->gcbsl1[p->level][q].cs>0)
#define QGCSL1LOOP QGCSLB1 QGCSLB1CHECK

#define QQGCSLB1 for(qq=0;qq<p->gcbsl1.ssize(p->level);++qq)
#define QQGCSLB1CHECK if(p->gcbsl1[p->level][qq].cs>0)
#define QQGCSL1LOOP QQGCSLB1 QQGCSLB1CHECK


#define GCSLB2 for(n=0;n<p->gcbsl2.ssize(p->level);++n)
#define GCSLB2CHECK if(p->gcbsl2[p->level][n].cs>0)
#define GCSL2LOOP GCSLB2 GCSLB2CHECK

#define QGCSLB2 for(q=0;q<p->gcbsl2.ssize(p->level);++q)
#define QGCSLB2CHECK if(p->gcbsl2[p->level][q].cs>0)
#define QGCSL2LOOP QGCSLB2 QGCSLB2CHECK

#define QQGCSLB2 for(qq=0;qq<p->gcbsl2.ssize(p->level);++qq)
#define QQGCSLB2CHECK if(p->gcbsl2[p->level][qq].cs>0)
#define QQGCSL2LOOP QQGCSLB2 QQGCSLB2CHECK


#define GCSLB4 for(n=0;n<p->gcbsl4.ssize(p->level);++n)
#define GCSLB4CHECK if(p->gcbsl4[p->level][n].cs>0)
#define GCSL4LOOP GCSLB4 GCSLB4CHECK

#define QGCSLB4 for(q=0;q<p->gcbsl4.ssize(p->level);++q)
#define QGCSLB4CHECK if(p->gcbsl4[p->level][q].cs>0)
#define QGCSL4LOOP QGCSLB4 QGCSLB4CHECK

#define QQGCSLB4 for(qq=0;qq<p->gcbsl4.ssize(p->level);++qq)
#define QQGCSLB4CHECK if(p->gcbsl4[p->level][qq].cs>0)
#define QQGCSL4LOOP QQGCSLB4 QQGCSLB4CHECK

// Inflow / outflow / active-wave-absorption lists, built in
// ghostcell::gcsl_setbcio. gcslin/gcslout come from gcbsl4, gcslawa1 from
// gcbsl1, gcslawa2 from gcbsl2.
#define GCSLIN for(n=0;n<p->gcslin.ssize(p->level);++n)
#define GCSLOUT for(n=0;n<p->gcslout.ssize(p->level);++n)
#define GCSLAWA1 for(n=0;n<p->gcslawa1.ssize(p->level);++n)
#define GCSLAWA2 for(n=0;n<p->gcslawa2.ssize(p->level);++n)


#define GCSLDFETA4 for(n=0;n<p->gcsldfeta4_count;++n)
#define GCSLDFETA4CHECK if(p->gcsldfeta4[n][2]>0)
#define GCSLDFETA4LOOP GCSLDFETA4 GCSLDFETA4CHECK

#define GCSLDFBED4 for(n=0;n<p->gcsldfbed4_count;++n)
#define GCSLDFBED4CHECK if(p->gcsldfbed4[n][2]>0)
#define GCSLDFBED4LOOP GCSLDFBED4 GCSLDFBED4CHECK

#endif
