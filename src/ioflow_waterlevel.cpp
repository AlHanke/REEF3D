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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void ioflow_f::fsfinflow(lexer *p, fdm *a, ghostcell *pgc)
{
    p->phimean=p->F56;

    if(p->F60>-1.0e20)
    {
        p->phimean=p->F60;
        p->phiout=p->F60;
    }

    if(p->F62>-1.0e20)
    p->phiout=p->F62;

    if(p->F61>-1.0e20)
    p->phiin=p->F62;

    count=0;
    zval=0.0;
    // LEVEL_LOOP
    GCINLOOP
    {
        i=p->gcin[p->level][n].i;
        j=p->gcin[p->level][n].j;
        k=p->gcin[p->level][n].k;

        GCIN_TILE(n);

        if(a->phi(i-1,j,k)>=0.0 && a->phi(i-1,j,k+1)<0.0)
        {
            zval+=-(a->phi(i-1,j,k)*p->DZP[KP])/(a->phi(i-1,j,k+1)-a->phi(i-1,j,k)) + p->pos_z();
            ++count;
        }
    }
    GC_TILE_RESET;

    count=pgc->globalisum(count);
    zval=pgc->globalsum(zval);

    if(count>0)
    {
        p->phimean=zval/double(count);

        if(p->F50==2 || p->F50==4)
        {
            // LEVEL_LOOP
            GCINLOOP
            {
                i=p->gcin[p->level][n].i;
                j=p->gcin[p->level][n].j;
                k=p->gcin[p->level][n].k;

                GCIN_TILE(n);

                a->phi(i-1,j,k)=p->phimean-p->pos_z();
                a->phi(i-2,j,k)=p->phimean-p->pos_z();
                a->phi(i-3,j,k)=p->phimean-p->pos_z();
            }
            GC_TILE_RESET;
        }
        p->phimean=pgc->globalmax(p->phimean);
    }

    // Outflow Water Level
    count=0;
    zval=0.0;
    // LEVEL_LOOP
    GCOUTLOOP
    {
        i=p->gcout[p->level][n].i;
        j=p->gcout[p->level][n].j;
        k=p->gcout[p->level][n].k;

        GCOUT_TILE(n);

        if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
        {
            zval+= -(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z();
            ++count;
        }
    }
    GC_TILE_RESET;

    count=pgc->globalisum(count);
    zval=pgc->globalsum(zval);

    if(count>0)
    {
        p->phiout=zval/double(count);
        p->phiout=pgc->globalmax(p->phiout);
    }

    // set outflow fsf
    double wsfout=p->phimean;
    double f;

    if(p->F62>1.0e-20)
    {
        if(p->F64==0)
        wsfout=p->F62;

        if(p->F64>0)
        {
            if(p->count<p->F64)
            f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;

            if(p->count>=p->F64)
            f = 1.0;

            wsfout = f*p->F62 + (1.0-f)*p->F60;
        }
    }

    if(p->F62>-1.0e20 && p->B77==2)
    {
        // LEVEL_LOOP
        GCOUTLOOP
        {
            i=p->gcout[p->level][n].i;
            j=p->gcout[p->level][n].j;
            k=p->gcout[p->level][n].k;

            GCOUT_TILE(n);

            a->phi(i+1,j,k)=wsfout-p->pos_z();
            a->phi(i+2,j,k)=wsfout-p->pos_z();
            a->phi(i+3,j,k)=wsfout-p->pos_z();
        }
        GC_TILE_RESET;
    }

    pBC->patchBC_waterlevel(p,a,pgc,a->phi);
}

void ioflow_f::fsfrkout(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    if(p->F62<-1.0e19 || p->B77!=2)
    {
        // LEVEL_LOOP
        GCOUTLOOP
        {
            i=p->gcout[p->level][n].i;
            j=p->gcout[p->level][n].j;
            k=p->gcout[p->level][n].k;

            GCOUT_TILE(n);

            f(i+1,j,k)=a->phi(i+1,j,k);
            f(i+2,j,k)=a->phi(i+2,j,k);
            f(i+3,j,k)=a->phi(i+3,j,k);
        }
        GC_TILE_RESET;
    }

    if(p->F62>-1.0e20 && p->B77==2)
    {
        // LEVEL_LOOP
        GCOUTLOOP
        {
            i=p->gcout[p->level][n].i;
            j=p->gcout[p->level][n].j;
            k=p->gcout[p->level][n].k;

            GCOUT_TILE(n);

            f(i+1,j,k)=p->F62-p->pos_z();
            f(i+2,j,k)=p->F62-p->pos_z();
            f(i+3,j,k)=p->F62-p->pos_z();
        }
        GC_TILE_RESET;
    }
}

void ioflow_f::fsfrkin(lexer *p, fdm *a, ghostcell *pgc, field& f)
{
    // LEVEL_LOOP
    GCINLOOP
    {
        i=p->gcin[p->level][n].i;
        j=p->gcin[p->level][n].j;
        k=p->gcin[p->level][n].k;

        GCIN_TILE(n);

        f(i-1,j,k)=a->phi(i-1,j,k);
        f(i-2,j,k)=a->phi(i-2,j,k);
        f(i-3,j,k)=a->phi(i-3,j,k);
    }
    GC_TILE_RESET;
}

double ioflow_f::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    return 0.0;
}

double ioflow_f::wave_xvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_f::wave_yvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

double ioflow_f::wave_zvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    return 0.0;
}

void ioflow_f::waterlevel_update(lexer *p,fdm *a,ghostcell *pgc)
{
    double zval;

    ILOOP
    JLOOP
    {
        zval=-1e20;

        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            zval=MAX(zval,-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
        }

        a->WL(i,j) = zval - a->bed(i,j);

        a->WL(i,j) = MAX(0.0001,a->WL(i,j));
    }
}
