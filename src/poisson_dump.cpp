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

#include "poisson_dump.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "field.h"
#include "looping.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <set>
#include <array>

void poisson_dump::start(lexer *p, fdm *a, ghostcell *pgc, field &pcorr, double alpha)
{
    const char *tag = std::getenv("REEF_POISSON_DUMP");
    if(!tag) return;

    // Stage counter: pjm_corr calls the solve once per RK stage, and the caller does not
    // pass the stage index. Count calls within a step instead, resetting whenever p->count
    // advances, so a dump can be pinned to one stage and compared like for like.
    static int last_count = -1;
    static int stage      = 0;
    if(p->count != last_count) { last_count = p->count; stage = 0; }
    else                       { ++stage; }

    if(const char *s = std::getenv("REEF_POISSON_DUMP_STEP"))
        if(p->count != std::atoi(s)) return;
    if(const char *s = std::getenv("REEF_POISSON_DUMP_STAGE"))
        if(stage != std::atoi(s)) return;

    // Metadata, rank 0, once. Records everything the analysis needs to know it is comparing
    // the same case: resolution per level, solver, tolerance and which build wrote this.
    if(p->mpirank == 0)
    {
        static bool meta_written = false;
        if(!meta_written)
        {
            meta_written = true;
            char mname[512];
            std::snprintf(mname, sizeof(mname), "pdump_%s_meta.txt", tag);
            std::ofstream m(mname);
            m << std::setprecision(17);
#if USE_AMREX
            m << "build USE_AMREX 1\n";
            m << "nlevs " << p->nlevs << "\n";
            for(int lev = 0; lev < p->nlevs; ++lev)
                m << "level " << lev
                  << " dx " << p->amrex_geometry[lev].CellSize(0)
                  << " dy " << p->amrex_geometry[lev].CellSize(1)
                  << " dz " << p->amrex_geometry[lev].CellSize(2) << "\n";
#else
            m << "build USE_AMREX 0\n";
            m << "nlevs 1\n";
            m << "level 0 dx " << p->dx << " dy " << p->dx << " dz " << p->dx << "\n";
#endif
            m << "N10 " << p->N10 << "\nN44 " << p->N44 << "\nN46 " << p->N46 << "\n";
            m << "W1 " << p->W1 << "\nW3 " << p->W3 << "\n";
            m << "psi " << p->psi << "\n";
            m << "DXM " << p->DXM << " DYM " << p->DYM << " DZM " << p->DZM << "\n";
            m << "F45 " << p->F45 << "\n";
            m << "j_dir " << p->j_dir << "\nmpi_ranks " << p->mpi_size << "\n";
            m << "cellnumtot " << p->cellnumtot << "\n";
            m.close();
        }
    }

#if USE_AMREX
    // Covered cells: amr_cell_mf is 1 on uncovered and 0 on covered for every level except
    // the finest, which is all 0 (see project_amr_cell_mf_convention -- a naive mask==0 test
    // inverts coverage there). Collect the covered set per level up front so the write loop
    // below stays a plain BASELOOP that compiles in both builds.
    std::set<std::array<int,4>> covered;
    for(int lev = 0; lev < p->nlevs - 1; ++lev)
    {
        const auto &cmf = p->amr_cell_mf[lev];
        for(amrex::MFIter mfi(cmf); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto ca = cmf.const_array(mfi);
            for(int kk = bx.smallEnd(2); kk <= bx.bigEnd(2); ++kk)
            for(int jj = bx.smallEnd(1); jj <= bx.bigEnd(1); ++jj)
            for(int ii = bx.smallEnd(0); ii <= bx.bigEnd(0); ++ii)
                if(ca(ii,jj,kk) == 0) covered.insert({lev,ii,jj,kk});
        }
    }
#endif

    char fname[512];
    std::snprintf(fname, sizeof(fname), "pdump_%s_s%05d_g%d_r%03d.dat",
                  tag, p->count, stage, p->mpirank);
    std::ofstream out(fname);
    out << std::setprecision(17) << std::scientific;
    out << "# lev i j k x y z flag4 phi ro rhs pcorr press covered\n";
    out << "# step " << p->count << " stage " << stage
        << " alpha " << alpha << " dt " << p->dt << " simtime " << p->simtime << "\n";

    BASELOOP
    {
        // Global index: BASELOOP's i,j,k are tile-local under AMReX. Emitted as metadata
        // only -- the join key is the coordinate triple, which is build-independent.
        int gi = i, gj = j, gk = k, lev = 0;
#if USE_AMREX
        gi += p->amr_tile_lo.x; gj += p->amr_tile_lo.y; gk += p->amr_tile_lo.z;
        lev = p->level;
        const int cov = covered.count({lev,gi,gj,gk}) ? 1 : 0;
#else
        const int cov = 0;
#endif
        // rhsvec is indexed by the row number poisson_pcorr stamped into Mrow, and is only
        // defined on solved (fluid) rows; everything else reports 0 rather than a stale slot.
        const int  row = a->Mrow(i,j,k);
        const double r = (row >= 0 && row < p->veclength) ? a->rhsvec.V[row] : 0.0;

        out << lev << ' ' << gi << ' ' << gj << ' ' << gk << ' '
            << p->pos_x() << ' ' << p->pos_y() << ' ' << p->pos_z() << ' '
            << p->flag4(i,j,k) << ' '
            << a->phi(i,j,k) << ' '
            << a->ro(i,j,k) << ' '
            << r << ' '
            << pcorr(i,j,k) << ' '
            << a->press(i,j,k) << ' '
            << cov << '\n';
    }
    out.close();
}
