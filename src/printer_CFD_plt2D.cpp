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

#include "printer_CFD.h"
#include "lexer.h"

#if USE_AMREX

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_FPC.H>
#include <AMReX_Utility.H>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

// Writes the pseudo-2D (ny==1) plot data as a genuine 2D (x,z) plotfile in the
// native AMReX format, readable by ParaView's AMReX/BoxLib reader with per-level
// selection. AMReX itself cannot write this from a 3D build (VisMF prints 3D boxes),
// but all dimension-dependent content is ASCII and the binary payload of an
// (nx,1,nz) FAB is byte-identical to an (nx,nz) one, so the three file types
// (Header, Cell_H, Cell_D) are emitted directly here.
void printer_CFD::print2D_plotfile_amrex(lexer* p,
                                         const amrex::Vector<amrex::MultiFab>& plot_mfs,
                                         const amrex::Vector<std::string>& varnames,
                                         int num)
{
    const int nlevs = p->nlevs;
    const int ncomp = (int)varnames.size();
    const int myproc = amrex::ParallelDescriptor::MyProc();
    const int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();

    const std::string pltname = amrex::Concatenate("REEF3D_CFD_PLT2D/plt", num, 7);

    if (amrex::ParallelDescriptor::IOProcessor())
        for (int lev = 0; lev < nlevs; lev++)
            amrex::UtilCreateDirectory(pltname + "/Level_" + std::to_string(lev), 0755);
    amrex::ParallelDescriptor::Barrier();

    // 2D box in native ASCII form, y dropped
    auto box2d = [](const amrex::Box& bx) {
        std::ostringstream os;
        os << "((" << bx.smallEnd(0) << ',' << bx.smallEnd(2) << ") ("
           << bx.bigEnd(0) << ',' << bx.bigEnd(2) << ") (0,0))";
        return os.str();
    };

    for (int lev = 0; lev < nlevs; lev++)
    {
        const amrex::MultiFab& mf = plot_mfs[lev];
        const amrex::BoxArray& ba = mf.boxArray();
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int nfabs = (int)ba.size();
        const std::string lvldir = pltname + "/Level_" + std::to_string(lev);

        // per-FAB metadata [byte offset, min[ncomp], max[ncomp]], owner fills its
        // slots, zeros elsewhere, so a sum-reduction collects everything on ioproc
        const int stride = 1 + 2*ncomp;
        std::vector<double> meta((size_t)nfabs*stride, 0.0);

        std::ofstream ofs;
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            if (!ofs.is_open())
                ofs.open(lvldir + "/" + amrex::Concatenate("Cell_D_", myproc, 5),
                         std::ios::out | std::ios::trunc | std::ios::binary);

            const amrex::Box& vbx = mfi.validbox();
            const int gid = mfi.index();
            const int ilo = vbx.smallEnd(0), ihi = vbx.bigEnd(0);
            const int klo = vbx.smallEnd(2), khi = vbx.bigEnd(2);
            const int jlo = vbx.smallEnd(1);
            const int nx = ihi - ilo + 1;

            meta[(size_t)gid*stride] = (double)(long long)ofs.tellp();

            ofs << "FAB " << amrex::FPC::NativeRealDescriptor()
                << box2d(vbx) << ' ' << ncomp << '\n';

            auto const& arr = mf.const_array(mfi);
            double* mn = &meta[(size_t)gid*stride + 1];
            double* mx = &meta[(size_t)gid*stride + 1 + ncomp];
            std::vector<double> row(nx);

            for (int c = 0; c < ncomp; c++)
            {
                mn[c] = std::numeric_limits<double>::max();
                mx[c] = std::numeric_limits<double>::lowest();
                for (int kk = klo; kk <= khi; kk++)
                {
                    for (int ii = ilo; ii <= ihi; ii++)
                    {
                        const double v = arr(ii,jlo,kk,c);
                        row[ii-ilo] = v;
                        mn[c] = std::min(mn[c], v);
                        mx[c] = std::max(mx[c], v);
                    }
                    ofs.write(reinterpret_cast<const char*>(row.data()), nx*sizeof(double));
                }
            }
        }
        if (ofs.is_open())
            ofs.close();

        amrex::ParallelDescriptor::ReduceRealSum(meta.data(), (int)meta.size(), ioproc);

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::ofstream hs(lvldir + "/Cell_H",
                             std::ios::out | std::ios::trunc | std::ios::binary);
            hs << 1 << '\n';        // VisMF::Header::Version_v1
            hs << 1 << '\n';        // How::NFiles
            hs << ncomp << '\n';
            hs << 0 << '\n';        // ngrow

            hs << '(' << nfabs << ' ' << 0 << '\n';
            for (int b = 0; b < nfabs; b++)
                hs << box2d(ba[b]) << '\n';
            hs << ")\n";

            hs << nfabs << '\n';
            for (int b = 0; b < nfabs; b++)
                hs << "FabOnDisk: " << amrex::Concatenate("Cell_D_", dm[b], 5)
                   << ' ' << (long long)meta[(size_t)b*stride] << '\n';
            hs << '\n';

            hs << std::scientific << std::setprecision(17);
            for (int m = 0; m < 2; m++) // min matrix, then max matrix
            {
                hs << nfabs << ',' << ncomp << '\n';
                for (int b = 0; b < nfabs; b++)
                {
                    for (int c = 0; c < ncomp; c++)
                        hs << meta[(size_t)b*stride + 1 + m*ncomp + c] << ',';
                    hs << '\n';
                }
                hs << '\n';
            }
        }
    }

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ofstream hs(pltname + "/Header",
                         std::ios::out | std::ios::trunc | std::ios::binary);
        hs.precision(17);

        hs << "HyperCLaw-V1.1\n";
        hs << ncomp << '\n';
        for (const auto& vn : varnames)
            hs << vn << '\n';
        hs << 2 << '\n';                    // spacedim
        hs << p->simtime << '\n';
        hs << nlevs-1 << '\n';

        // problem extents in output coordinates, consistent with the 3D plotfile path
        const auto& rb0 = p->amrex_geometry[0].ProbDomain();
        const double xlo_out = p->coordinates::Xout(rb0.lo(0), rb0.lo(1));
        const double xhi_out = p->coordinates::Xout(rb0.hi(0), rb0.hi(1));
        hs << xlo_out << ' ' << rb0.lo(2) << " \n";
        hs << xhi_out << ' ' << rb0.hi(2) << " \n";

        for (int lev = 0; lev < nlevs-1; lev++)
            hs << p->ref_vec[0] << ' ';
        hs << '\n';
        for (int lev = 0; lev < nlevs; lev++)
            hs << box2d(p->amrex_geometry[lev].Domain()) << ' ';
        hs << '\n';
        for (int lev = 0; lev < nlevs; lev++)
            hs << p->count << ' ';
        hs << '\n';
        for (int lev = 0; lev < nlevs; lev++)
            hs << p->amrex_geometry[lev].CellSize(0) << ' '
               << p->amrex_geometry[lev].CellSize(2) << " \n";
        hs << 0 << '\n';                    // CoordSys::cartesian
        hs << 0 << '\n';                    // boundary width

        for (int lev = 0; lev < nlevs; lev++)
        {
            const amrex::BoxArray& ba = plot_mfs[lev].boxArray();
            const double dx = p->amrex_geometry[lev].CellSize(0);
            const double dz = p->amrex_geometry[lev].CellSize(2);

            hs << lev << ' ' << ba.size() << ' ' << p->simtime << '\n';
            hs << p->count << '\n';
            for (int b = 0; b < (int)ba.size(); b++)
            {
                const amrex::Box& bx = ba[b];
                hs << xlo_out + bx.smallEnd(0)*dx << ' '
                   << xlo_out + (bx.bigEnd(0)+1)*dx << '\n';
                hs << rb0.lo(2) + bx.smallEnd(2)*dz << ' '
                   << rb0.lo(2) + (bx.bigEnd(2)+1)*dz << '\n';
            }
            hs << "Level_" << lev << "/Cell\n";
        }
    }

    amrex::ParallelDescriptor::Barrier();
}

#endif
