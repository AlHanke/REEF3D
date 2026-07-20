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

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::gcsl_setbc1(lexer *p)
{
    int cs,bc;

    GCSL1LOOP
    {
        i = p->gcbsl1[p->level][n].i;
        j = p->gcbsl1[p->level][n].j;
        cs = p->gcbsl1[p->level][n].cs;
        bc = p->gcbsl1[p->level][n].bc;

        if(bc==WALL)
        {
            if(cs==X_NEG && i+p->origin_i==0)
                p->gcbsl1[p->level][n].bc=p->bcside1;
            else if(cs==X_POS && i+p->origin_i==p->gknox-2)
                p->gcbsl1[p->level][n].bc=p->bcside4;
            else if(cs==Y_NEG && j+p->origin_j==0)
                p->gcbsl1[p->level][n].bc=p->bcside3;
            else if(cs==Y_POS && j+p->origin_j==p->gknoy-1)
                p->gcbsl1[p->level][n].bc=p->bcside2;
        }
    }
}

void ghostcell::gcsl_setbc2(lexer *p)
{
    int cs,bc;

    GCSL2LOOP
    {
        i = p->gcbsl2[p->level][n].i;
        j = p->gcbsl2[p->level][n].j;
        cs = p->gcbsl2[p->level][n].cs;
        bc = p->gcbsl2[p->level][n].bc;

        if(bc==WALL)
        {
            if(cs==X_NEG && i+p->origin_i==0)
                p->gcbsl2[p->level][n].bc=p->bcside1;
            else if(cs==X_POS && i+p->origin_i==p->gknox-1)
                p->gcbsl2[p->level][n].bc=p->bcside4;
            else if(cs==Y_NEG && j+p->origin_j==0)
                p->gcbsl2[p->level][n].bc=p->bcside3;
            else if(cs==Y_POS && j+p->origin_j==p->gknoy-2)
                p->gcbsl2[p->level][n].bc=p->bcside2;
        }
    }
}

void ghostcell::gcsl_setbc4(lexer *p)
{
    int cs,bc;

    GCSL4LOOP
    {
        i = p->gcbsl4[p->level][n].i;
        j = p->gcbsl4[p->level][n].j;
        cs = p->gcbsl4[p->level][n].cs;
        bc = p->gcbsl4[p->level][n].bc;

        if(bc==WALL)
        {
            if(cs==X_NEG && i+p->origin_i==0)
                p->gcbsl4[p->level][n].bc=p->bcside1;
            else if(cs==X_POS && i+p->origin_i==p->gknox-1)
                p->gcbsl4[p->level][n].bc=p->bcside4;
            else if(cs==Y_NEG && j+p->origin_j==0)
                p->gcbsl4[p->level][n].bc=p->bcside3;
            else if(cs==Y_POS && j+p->origin_j==p->gknoy-1)
                p->gcbsl4[p->level][n].bc=p->bcside2;
        }
    }
}

void ghostcell::gcsl_setbcio(lexer *p)
{
    int cs,bc;

    // Each list sizes itself from its own source list. gcslawa1/gcslawa2 are
    // built from gcbsl1/gcbsl2, so they must not be sized off gcbsl4's
    // outflow count.
    p->gcslin.resize_levels(p->gcbsl4.nlevels());
    p->gcslout.resize_levels(p->gcbsl4.nlevels());

    // gcslawa1/gcslawa2 is not used by FNPF so gcbsl1 & gcbsl2 not being build
    // is not a problem.
    p->gcslawa1.resize_levels(p->gcbsl1.nlevels());
    p->gcslawa2.resize_levels(p->gcbsl2.nlevels());

    p->gcslin[p->level].clear();
    p->gcslout[p->level].clear();

    GCSL4LOOP
    {
        i = p->gcbsl4[p->level][n].i;
        j = p->gcbsl4[p->level][n].j;
        cs = p->gcbsl4[p->level][n].cs;
        bc = p->gcbsl4[p->level][n].bc;

        if(bc==INFLOW || bc==WAVEGEN)
        {
            auto &e = p->gcslin[p->level].emplace_back(gcb_sl{i,j});
            GCB_COPY_TILE(e, p->gcbsl4[p->level][n]);
        }
        else if(bc==OUTFLOW || bc==NUMBEACH)
        {
            auto &e = p->gcslout[p->level].emplace_back(gcb_sl_cs{i,j,cs});
            GCB_COPY_TILE(e, p->gcbsl4[p->level][n]);
        }
    }

    p->gcslawa1[p->level].clear();
    GCSL1LOOP
    {
        i = p->gcbsl1[p->level][n].i;
        j = p->gcbsl1[p->level][n].j;
        cs = p->gcbsl1[p->level][n].cs;
        bc = p->gcbsl1[p->level][n].bc;

        if(bc==OUTFLOW || bc==NUMBEACH)
        {
            auto &e = p->gcslawa1[p->level].emplace_back(gcb_sl_cs{i,j,cs});
            GCB_COPY_TILE(e, p->gcbsl1[p->level][n]);
        }
    }

    p->gcslawa2[p->level].clear();
    GCSL2LOOP
    {
        i = p->gcbsl2[p->level][n].i;
        j = p->gcbsl2[p->level][n].j;
        cs = p->gcbsl2[p->level][n].cs;
        bc = p->gcbsl2[p->level][n].bc;

        if(bc==OUTFLOW || bc==NUMBEACH)
        {
            auto &e = p->gcslawa2[p->level].emplace_back(gcb_sl_cs{i,j,cs});
            GCB_COPY_TILE(e, p->gcbsl2[p->level][n]);
        }
    }

    // IOSL
    p->IOSL.setVal(0, true);

    // IOSL is an ArrayWrapper2D: operator() resolves through the
    // installed tile context, so unlike the derivation loops above this one has
    // to reinstate the tile each entry was recorded under. Safe against
    // GCSL4LOOP's bound (which re-reads p->level, and set_tile_ctx writes it)
    // because the entry's context level is the level its list is stored under.
    GCSL4LOOP
    {
        GCB_TILE(p->gcbsl4[p->level][n], p->level);

        i = p->gcbsl4[p->level][n].i;
        j = p->gcbsl4[p->level][n].j;
        cs = p->gcbsl4[p->level][n].cs;
        bc = p->gcbsl4[p->level][n].bc;

        if((bc==INFLOW || bc==WAVEGEN) && cs==X_NEG)
            p->IOSL(i-1,j)=1;
        else if((bc==OUTFLOW || bc==NUMBEACH) && cs==X_POS)
            p->IOSL(i+1,j)=2;
    }
    GC_TILE_RESET;
}
