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

#include "poisson_pcorr.h"
#include "lexer.h"
#include "fdm.h"
#include "heat.h"
#include "concentration.h"
#include "density_f.h"
#include "density_sf.h"
#include "density_comp.h"
#include "density_conc.h"
#include "density_heat.h"
#include "density_vof.h"
#include "density_rheo.h"
#include "density_pst.h"
#include <cstdlib>

poisson_pcorr::poisson_pcorr(lexer *p, heat *&pheat, concentration *&pconc)
{
    if(p->F80==0 && p->F300==0 && p->W90==0)
    {
        if(p->W30==0 && p->C10==0 && p->H10==0)
        {
            if(p->Q10==0)
            pd = new density_f(p);
            else if(p->Q10==1)
            pd = new density_f(p);
        }

        if(p->H10==0 && p->W30==1)
        pd = new density_comp(p);

        if(p->H10>0 && p->C10==0)
        pd = new density_heat(p,pheat);

        if(p->C10>0 && p->H10==0)
        pd = new density_conc(p,pconc);
    }

    if(p->F80>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0)
    pd = new density_vof(p);

    if(p->F30>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);

    if(p->F300>=1)
    pd = new density_rheo(p);
}

poisson_pcorr::~poisson_pcorr()
{
    delete pd;
}

void poisson_pcorr::start(lexer* p, fdm *a, field &press)
{
    a->M.reset();

    const bool is3D = p->j_dir;

    // Materialise rho_face into a->rofx/rofy/rofz once, instead of evaluating the density
    // model six times per cell below. With AMR levels this is also where the C-F faces are
    // made single-valued (density::update_faces -> amrex::average_down_faces), so the coarse
    // and fine sides of one physical face can no longer disagree when the smoothed Heaviside
    // band crosses the interface. Cell (i,j,k) stores its HIGH face, so a low face reads
    // (i-1,j,k).
    // TODO: hoist to a single per-RK-stage call once the predictor (pjm_corr u/v/wpgrad) and
    // hypre_ssamg::velcorr_amr read the same arrays; it lives here while poisson_pcorr is the
    // only converted consumer.
    pd->update_faces(p,a);

    // Row numbering: stamp a->Mrow with the LOOP index used to fill a->M / a->rhsvec.
    // This is the single source of truth consumed by hypre_ssamg (matrix fill) and the
    // velocity corrector -- keep it on the coefficient LOOP so index n and the a->M[n]
    // writes below are defined in one pass and can never drift.
    n=0;
    LOOP
    {
        a->Mrow(i,j,k) = n;

        a->M.p[n] = (CPOR1*PORVAL1)/(a->rofx(i,j,k)*p->DXP[IP]*p->DXN[IP])
                  + (CPOR1m*PORVAL1m)/(a->rofx(i-1,j,k)*p->DXP[IM1]*p->DXN[IP])

                  +(is3D ? ((CPOR2*PORVAL2)/(a->rofy(i,j,k)*p->DYP[JP]*p->DYN[JP])
                  + (CPOR2m*PORVAL2m)/(a->rofy(i,j-1,k)*p->DYP[JM1]*p->DYN[JP])) : 0.0)

                  + (CPOR3*PORVAL3)/(a->rofz(i,j,k)*p->DZP[KP]*p->DZN[KP])
                  + (CPOR3m*PORVAL3m)/(a->rofz(i,j,k-1)*p->DZP[KM1]*p->DZN[KP]);


        a->M.n[n] = -(CPOR1*PORVAL1)/(a->rofx(i,j,k)*p->DXP[IP]*p->DXN[IP]);
        a->M.s[n] = -(CPOR1m*PORVAL1m)/(a->rofx(i-1,j,k)*p->DXP[IM1]*p->DXN[IP]);

        if(is3D)
        {
            a->M.w[n] = -(CPOR2*PORVAL2)/(a->rofy(i,j,k)*p->DYP[JP]*p->DYN[JP]);
            a->M.e[n] = -(CPOR2m*PORVAL2m)/(a->rofy(i,j-1,k)*p->DYP[JM1]*p->DYN[JP]);
        }
        else
        {
            a->M.w[n] = 0.0;
            a->M.e[n] = 0.0;
        }

        a->M.t[n] = -(CPOR3*PORVAL3)/(a->rofz(i,j,k)*p->DZP[KP]*p->DZN[KP]);
        a->M.b[n] = -(CPOR3m*PORVAL3m)/(a->rofz(i,j,k-1)*p->DZP[KM1]*p->DZN[KP]);

        ++n;
    }

    const bool wallfold = (std::getenv("REEF_NO_WALLFOLD")==nullptr); // TEST: set env to disable solidnb fold
#if USE_AMREX
    const bool symfold = (std::getenv("REEF_SYM_FOLD")!=nullptr);
#endif
    LOOP
    {
        const int n = a->Mrow(i,j,k);
        if(n < 0) continue;

        // ---- Solid-wall faces (Option B): fold into the diagonal (homogeneous Neumann) every
        // face whose neighbour is a SOLID. This ties the matrix coupling structure to the velocity
        // DOFs (the correction skips non-DOF faces), giving D(1/rho)G = L at walls and leaving every
        // retained coupling a real velocity face with two valid pressure cells. Folding zeroes
        // M.face, so the open-boundary Dirichlet/outflow branches below no-op on already-folded faces.
        //
        // "Solid" is read from the NEIGHBOUR cell flag4, NOT the face flag1/2/3: a low face carries
        // the neighbour's flag (flag1(i-1)=flag4(i-1)), but a HIGH face is forced to -10 regardless
        // of neighbour type, so flag1/2/3 cannot tell a solid wall from the free surface. flag4 can.
        //
        // EXCLUDED, kept Dirichlet: AIR_FLAG -- the free surface is also a non-DOF face, but its
        // atmospheric Dirichlet pressure is the ONLY anchor of this closed domain; folding it leaves
        // an all-Neumann (singular) operator and the solve diverges. Also INFLOW/OUTFLOW/IO==2.
        // (Domain solid walls currently carry the AIR_FLAG ghost too, so they stay Dirichlet here;
        //  giving them a real solid flag is the remaining step to fold them as well.)
        auto solidnb = [&](int f, int io){ return wallfold && f<0 && f!=AIR_FLAG && f!=INFLOW_FLAG && f!=OUTFLOW_FLAG && io!=2; };

        if(solidnb(p->flag4(i-1,j,k), p->IO(i-1,j,k))) { a->M.p[n] += a->M.s[n]; a->M.s[n] = 0.0; }
        if(solidnb(p->flag4(i+1,j,k), p->IO(i+1,j,k))) { a->M.p[n] += a->M.n[n]; a->M.n[n] = 0.0; }

        if(is3D)
        {
            if(solidnb(p->flag4(i,j-1,k), p->IO(i,j-1,k))) { a->M.p[n] += a->M.e[n]; a->M.e[n] = 0.0; }
            if(solidnb(p->flag4(i,j+1,k), p->IO(i,j+1,k))) { a->M.p[n] += a->M.w[n]; a->M.w[n] = 0.0; }
        }

        if(solidnb(p->flag4(i,j,k-1), p->IO(i,j,k-1))) { a->M.p[n] += a->M.b[n]; a->M.b[n] = 0.0; }
        if(solidnb(p->flag4(i,j,k+1), p->IO(i,j,k+1))) { a->M.p[n] += a->M.t[n]; a->M.t[n] = 0.0; }

#if USE_AMREX
        // ---- TEST (env REEF_SYM_FOLD): geometry-based fold. Fold every face that sits on the
        // level's DOMAIN BOUNDARY into the diagonal (homogeneous Neumann), so symmetry planes get
        // the same pcorr BC as solid walls WITHOUT relabeling flag4. The free surface is interior
        // (never on the domain boundary) so it is untouched and keeps its Dirichlet anchor.
        // Periodic boundaries and outflow are excluded. Already-folded faces (solid bottom) no-op
        // since M.face is 0. Global index = local + amr_tile_lo, compared to the level's Domain().
        if(symfold)
        {
            const auto dom = p->amrex_geometry[p->level].Domain();
            const int gi=i+p->amr_tile_lo.x, gj=j+p->amr_tile_lo.y, gk=k+p->amr_tile_lo.z;

            if(gi==dom.smallEnd(0) && p->periodic1==0 && p->flag4(i-1,j,k)!=OUTFLOW_FLAG && p->IO(i-1,j,k)!=2)
            { a->M.p[n] += a->M.s[n]; a->M.s[n] = 0.0; }
            if(gi==dom.bigEnd(0)   && p->periodic1==0 && p->flag4(i+1,j,k)!=OUTFLOW_FLAG && p->IO(i+1,j,k)!=2)
            { a->M.p[n] += a->M.n[n]; a->M.n[n] = 0.0; }

            if(is3D)
            {
                if(gj==dom.smallEnd(1) && p->periodic2==0 && p->flag4(i,j-1,k)!=OUTFLOW_FLAG && p->IO(i,j-1,k)!=2)
                { a->M.p[n] += a->M.e[n]; a->M.e[n] = 0.0; }
                if(gj==dom.bigEnd(1)   && p->periodic2==0 && p->flag4(i,j+1,k)!=OUTFLOW_FLAG && p->IO(i,j+1,k)!=2)
                { a->M.p[n] += a->M.w[n]; a->M.w[n] = 0.0; }
            }

            if(gk==dom.smallEnd(2) && p->periodic3==0 && p->flag4(i,j,k-1)!=OUTFLOW_FLAG && p->IO(i,j,k-1)!=2)
            { a->M.p[n] += a->M.b[n]; a->M.b[n] = 0.0; }
            if(gk==dom.bigEnd(2)   && p->periodic3==0 && p->flag4(i,j,k+1)!=OUTFLOW_FLAG && p->IO(i,j,k+1)!=2)
            { a->M.p[n] += a->M.t[n]; a->M.t[n] = 0.0; }
        }
#endif

        // ---- x direction
        // inflow
        if(p->flag4(i-1,j,k)<0 && (i+p->origin_i>0 || p->periodic1==0))
        {
            a->rhsvec.V[n] -= a->M.s[n]*press(i-1,j,k);
            a->M.s[n] = 0.0;
        }

        // outflow
        if(p->flag4(i+1,j,k)<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && (p->IO(i+1,j,k)!=2 || (p->B60!=1&&p->B99<3)))
        {
            a->rhsvec.V[n] -= a->M.n[n]*press(i+1,j,k);
            a->M.n[n] = 0.0;
        }

        // controlled outflow
        if(p->IO(i+1,j,k)==2)
        {
            if(p->B77==1)
            {
                if(p->F50==2 || p->F50==3)
                pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
                else if(p->F50==1 || p->F50==4)
                pval=a->press(i,j,k);
            }
            else if(p->B77==2)
            {
                pval=a->press(i,j,k);
            }
            else if(p->B77==10)
            pval=0.0;

            a->rhsvec.V[n] -= a->M.n[n]*(-a->press(i,j,k)+pval);
            a->M.n[n] = 0.0;
        }

        // AWA outflow
        if(p->flag4(i+1,j,k)<0 && (i+p->origin_i<p->gknox-1 || p->periodic1==0) && (p->IO(i+1,j,k)==2 && p->B90==1 && p->B99>2))
        {
            pval = (p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

            a->rhsvec.V[n] -= a->M.n[n]*(-a->press(i,j,k)+pval);
            a->M.n[n] = 0.0;
        }

        // ---- y direction
        if(is3D)
        {
            // east face
            if(p->flag4(i,j-1,k)<0 && (j+p->origin_j>0 || p->periodic2==0))
            {
                a->rhsvec.V[n] -= a->M.e[n]*press(i,j-1,k);
                a->M.e[n] = 0.0;
            }
            // west face
            if(p->flag4(i,j+1,k)<0 && (j+p->origin_j<p->gknoy-1 || p->periodic2==0))
            {
                a->rhsvec.V[n] -= a->M.w[n]*press(i,j+1,k);
                a->M.w[n] = 0.0;
            }
        }

        // ---- z direction
        // bottom face
        if(p->flag4(i,j,k-1)<0 && (k+p->origin_k>0 || p->periodic3==0))
        {
            a->rhsvec.V[n] -= a->M.b[n]*press(i,j,k-1);
            a->M.b[n] = 0.0;
        }
        // top face
        if(p->flag4(i,j,k+1)<0 && (k+p->origin_k<p->gknoz-1 || p->periodic3==0))
        {
            a->rhsvec.V[n] -= a->M.t[n]*press(i,j,k+1);
            a->M.t[n] = 0.0;
        }
    }
}
