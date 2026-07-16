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
#include "fdm.h"

#if USE_AMREX

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
    // one raw-appended binary block in a .vtu file: uint64 byte count + payload
    template<typename T>
    void append_block(std::ofstream& f, const std::vector<T>& v)
    {
        const uint64_t nbytes = v.size()*sizeof(T);
        f.write(reinterpret_cast<const char*>(&nbytes), sizeof(uint64_t));
        f.write(reinterpret_cast<const char*>(v.data()), nbytes);
    }

    // maintain a ParaView .pvd collection. The file is rewritten each print with
    // stale entries dropped — an entry is stale if it references the same dataset
    // file or a timestep at/after the new one — so restarted or repeated runs
    // never accumulate duplicate DataSet entries.
    void append_pvd(const std::string& pvd_path, const std::string& dataset, double time)
    {
        std::vector<std::string> kept;
        std::ifstream in(pvd_path);
        std::string line;
        while (std::getline(in, line))
        {
            const auto tpos = line.find("timestep=\"");
            const auto fpos = line.find("file=\"");
            if (tpos == std::string::npos || fpos == std::string::npos)
                continue; // header/footer line
            const double t = atof(line.c_str() + tpos + 10);
            const auto fend = line.find('"', fpos + 6);
            const std::string f = line.substr(fpos + 6, fend - fpos - 6);
            if (t < time && f != dataset)
                kept.push_back(line);
        }
        in.close();

        std::ofstream pvd(pvd_path, std::ios::out | std::ios::trunc | std::ios::binary);
        pvd << "<?xml version=\"1.0\"?>\n"
               "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
               "  <Collection>\n";
        for (const auto& l : kept)
            pvd << l << '\n';
        char entry[400];
        snprintf(entry, sizeof(entry),
                 "    <DataSet timestep=\"%.17g\" part=\"0\" file=\"%s\"/>\n",
                 time, dataset.c_str());
        pvd << entry
            << "  </Collection>\n"
               "</VTKFile>\n";
        pvd.close();
    }
}

// Point-interpolated output with AMR level structure as a vtkMultiBlockDataSet:
// one binary .vtu per (level, rank) holding the level's uncovered cells with
// node-interpolated point data, grouped into per-level blocks of a .vtm
// collection so ParaView's block selection acts as level selection; all blocks
// together form the composite (finest-available) field. Velocities are
// interpolated straight from the staggered faces (REEF3D convention:
// u(i,j,k) = east face of cell (i,j,k), v north, w top), scalars from the
// ghosted cell data, matching the legacy pre-AMReX point output. plane2D
// collapses the pseudo-2D y direction to a quad mesh in the y=0 plane.
// Alternative to print_interp_amrex (vtkPartitionedDataSetCollection): .vtm is
// the legacy composite format, but its reader enumerates data arrays for
// read-time selection, which ParaView's PDC reader does not (as of 6.1.1).
void printer_CFD::print_interp_amrex_vtm(lexer* p, fdm* a,
                                         const amrex::Vector<amrex::MultiFab>& plot_mfs,
                                         const amrex::Vector<std::string>& varnames,
                                         int num, bool plane2D)
{
    const int nlevs = p->nlevs;
    const int ncomp = (int)varnames.size();
    const int myproc = amrex::ParallelDescriptor::MyProc();
    const int nprocs = amrex::ParallelDescriptor::NProcs();

    const std::string base = plane2D ? "REEF3D_CFD_VTM2D" : "REEF3D_CFD_VTM";
    char stepname[100];
    snprintf(stepname, sizeof(stepname), "plt%07d", num);

    if (amrex::ParallelDescriptor::IOProcessor())
        amrex::UtilCreateDirectory(base + "/" + stepname, 0755);
    amrex::ParallelDescriptor::Barrier();

    int elev_comp = -1; // "elevation" is exact at nodes, no averaging
    for (int c = 0; c < ncomp; c++)
        if (varnames[c] == "elevation")
            elev_comp = c;

    for (int lev = 0; lev < nlevs; lev++)
    {
        const amrex::MultiFab& data_mf = plot_mfs[lev];
        const amrex::iMultiFab& mask_mf = p->amr_cell_mf[lev];
        const amrex::MultiFab& u_mf = a->u.GetMultiFab(lev);
        const amrex::MultiFab& v_mf = a->v.GetMultiFab(lev);
        const amrex::MultiFab& w_mf = a->w.GetMultiFab(lev);

        const amrex::Box dom = p->amrex_geometry[lev].Domain();
        const double dx = p->amrex_geometry[lev].CellSize(0);
        const double dy = p->amrex_geometry[lev].CellSize(1);
        const double dz = p->amrex_geometry[lev].CellSize(2);
        const double x0 = p->amrex_geometry[lev].ProbLo(0);
        const double y0 = p->amrex_geometry[lev].ProbLo(1);
        const double z0 = p->amrex_geometry[lev].ProbLo(2);

        // amr_cell_mf convention: the finest level is all 0; coarser levels are
        // built with makeFineMask(...,1,0) -> uncovered = 1, covered by fine = 0
        const bool finest = (lev == nlevs-1);

        // mesh + point data accumulated over all local boxes of this level
        std::unordered_map<uint64_t,int32_t> nodeid;
        std::vector<double> pts;
        std::vector<float> vel;
        std::vector<std::vector<float>> vals(ncomp);
        std::vector<int64_t> conn, offs;
        std::vector<uint8_t> types;

        for (amrex::MFIter mfi(data_mf); mfi.isValid(); ++mfi)
        {
            const amrex::Box& vbx = mfi.validbox();
            auto const& mask = mask_mf.const_array(mfi);
            auto const& cc = data_mf.const_array(mfi);
            auto const& uf = u_mf.const_array(mfi);
            auto const& vf = v_mf.const_array(mfi);
            auto const& wf = w_mf.const_array(mfi);

            // fetch with corner-ghost protection: values outside the level domain
            // in one direction are BC/C-F-filled ghosts and used as-is (legacy ipol
            // behavior); outside in >=2 directions (corner/edge ghosts, not reliably
            // filled) the offending indices are clamped, i.e. mirrored onto the domain.
            // In pseudo-2D runs the solver never fills y-direction ghosts at all, so
            // the y index is always collapsed onto the domain (legacy ipol forced j=0).
            const bool thin_y = !p->j_dir;
            auto out_lo = [&](int c, int d){ return c < dom.smallEnd(d); };
            auto out_hi = [&](int c, int d){ return c > dom.bigEnd(d); };
            auto nouts = [&](int ii, int jj, int kk){
                return int(out_lo(ii,0)||out_hi(ii,0)) + int(out_lo(jj,1)||out_hi(jj,1))
                     + int(out_lo(kk,2)||out_hi(kk,2)); };
            auto clampd = [&](int c, int d){
                return std::min(std::max(c, dom.smallEnd(d)), dom.bigEnd(d)); };

            auto cell = [&](int ii, int jj, int kk, int c){
                if (thin_y) jj = clampd(jj,1);
                if (nouts(ii,jj,kk) >= 2)
                    { ii = clampd(ii,0); jj = clampd(jj,1); kk = clampd(kk,2); }
                return cc(ii,jj,kk,c); };
            // staggered fetch: 'norm' is the face-normal direction, where index
            // dom.smallEnd(norm)-1 is the physical boundary face, not a corner ghost
            auto face = [&](const amrex::Array4<const amrex::Real>& f,
                            int ii, int jj, int kk, int norm){
                if (thin_y) jj = clampd(jj,1);
                const int idx[3] = {ii,jj,kk};
                int outs = 0;
                for (int d = 0; d < 3; d++)
                {
                    if (d == norm) outs += int(idx[d] < dom.smallEnd(d)-1 || idx[d] > dom.bigEnd(d));
                    else           outs += int(out_lo(idx[d],d) || out_hi(idx[d],d));
                }
                if (outs >= 2)
                {
                    int ci = ii, cj = jj, ck = kk;
                    if (norm != 0) ci = clampd(ii,0);
                    if (norm != 1) cj = clampd(jj,1);
                    if (norm != 2) ck = clampd(kk,2);
                    return f(ci,cj,ck);
                }
                return f(ii,jj,kk); };

            // create (or look up) the node at corner (i,j,k) and compute its values
            auto getnode = [&](int i, int j, int k) -> int32_t
            {
                const uint64_t key = (uint64_t(uint32_t(i)) << 42)
                                   | (uint64_t(uint32_t(j)) << 21)
                                   |  uint64_t(uint32_t(k));
                auto it = nodeid.find(key);
                if (it != nodeid.end())
                    return it->second;

                const int32_t id = (int32_t)(pts.size()/3);
                nodeid.emplace(key, id);

                const double px = x0 + i*dx, py = y0 + j*dy, pz = z0 + k*dz;
                pts.push_back(p->coordinates::Xout(px,py));
                pts.push_back(plane2D ? 0.0 : p->coordinates::Yout(px,py));
                pts.push_back(pz);

                double un, vn, wn;
                if (plane2D)
                {
                    un = 0.5 *( face(uf,i-1,j  ,k-1,0) + face(uf,i-1,j  ,k  ,0) );
                    vn = 0.25*( face(vf,i-1,j-1,k-1,1) + face(vf,i  ,j-1,k-1,1)
                              + face(vf,i-1,j-1,k  ,1) + face(vf,i  ,j-1,k  ,1) );
                    wn = 0.5 *( face(wf,i-1,j  ,k-1,2) + face(wf,i  ,j  ,k-1,2) );
                }
                else
                {
                    un = 0.25*( face(uf,i-1,j-1,k-1,0) + face(uf,i-1,j  ,k-1,0)
                              + face(uf,i-1,j-1,k  ,0) + face(uf,i-1,j  ,k  ,0) );
                    vn = 0.25*( face(vf,i-1,j-1,k-1,1) + face(vf,i  ,j-1,k-1,1)
                              + face(vf,i-1,j-1,k  ,1) + face(vf,i  ,j-1,k  ,1) );
                    wn = 0.25*( face(wf,i-1,j-1,k-1,2) + face(wf,i  ,j-1,k-1,2)
                              + face(wf,i-1,j  ,k-1,2) + face(wf,i  ,j  ,k-1,2) );
                }
                vel.push_back((float)un);
                vel.push_back((float)vn);
                vel.push_back((float)wn);
                vals[0].push_back((float)un);
                vals[1].push_back((float)vn);
                vals[2].push_back((float)wn);

                for (int c = 3; c < ncomp; c++)
                {
                    double s;
                    if (c == elev_comp)
                        s = pz;
                    else if (plane2D)
                        s = 0.25*( cell(i-1,j,k-1,c) + cell(i,j,k-1,c)
                                 + cell(i-1,j,k  ,c) + cell(i,j,k  ,c) );
                    else
                        s = 0.125*( cell(i-1,j-1,k-1,c) + cell(i,j-1,k-1,c)
                                  + cell(i-1,j  ,k-1,c) + cell(i,j  ,k-1,c)
                                  + cell(i-1,j-1,k  ,c) + cell(i,j-1,k  ,c)
                                  + cell(i-1,j  ,k  ,c) + cell(i,j  ,k  ,c) );
                    vals[c].push_back((float)s);
                }
                return id;
            };

            const int jplane = vbx.smallEnd(1);

            amrex::Loop(vbx, [&](int i, int j, int k)
            {
                if (!finest && mask(i,j,k) == 0) return;   // covered by a finer level
                if (plane2D && j != jplane) return;

                if (plane2D)
                {
                    conn.push_back(getnode(i  ,j,k  ));
                    conn.push_back(getnode(i+1,j,k  ));
                    conn.push_back(getnode(i+1,j,k+1));
                    conn.push_back(getnode(i  ,j,k+1));
                    types.push_back(9);            // VTK_QUAD
                }
                else
                {
                    conn.push_back(getnode(i  ,j  ,k  ));
                    conn.push_back(getnode(i+1,j  ,k  ));
                    conn.push_back(getnode(i+1,j+1,k  ));
                    conn.push_back(getnode(i  ,j+1,k  ));
                    conn.push_back(getnode(i  ,j  ,k+1));
                    conn.push_back(getnode(i+1,j  ,k+1));
                    conn.push_back(getnode(i+1,j+1,k+1));
                    conn.push_back(getnode(i  ,j+1,k+1));
                    types.push_back(12);           // VTK_HEXAHEDRON
                }
                offs.push_back((int64_t)conn.size());
            });
        }

        // ---- write this (level, rank) piece; empty pieces are valid VTK and keep
        //      the .vtm index deterministic across ranks
        const int64_t npts = (int64_t)(pts.size()/3);
        const int64_t ncells = (int64_t)types.size();

        char vtu_name[300];
        snprintf(vtu_name, sizeof(vtu_name), "%s/%s/lev%d_r%04d.vtu",
                 base.c_str(), stepname, lev, myproc);

        std::ofstream f(vtu_name, std::ios::out | std::ios::trunc | std::ios::binary);
        f << "<?xml version=\"1.0\"?>\n"
             "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
             "  <UnstructuredGrid>\n"
             "    <FieldData>\n"
             "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\"> "
          << std::scientific << p->simtime << " </DataArray>\n"
             "    </FieldData>\n"
             "    <Piece NumberOfPoints=\"" << npts << "\" NumberOfCells=\"" << ncells << "\">\n";

        uint64_t off = 0;
        auto blocksize = [](uint64_t nbytes){ return sizeof(uint64_t) + nbytes; };

        f << "      <Points>\n"
             "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << off << "\"/>\n"
             "      </Points>\n";
        off += blocksize(pts.size()*sizeof(double));

        f << "      <Cells>\n"
             "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"" << off << "\"/>\n";
        off += blocksize(conn.size()*sizeof(int64_t));
        f << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"" << off << "\"/>\n";
        off += blocksize(offs.size()*sizeof(int64_t));
        f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << off << "\"/>\n";
        off += blocksize(types.size()*sizeof(uint8_t));
        f << "      </Cells>\n";

        f << "      <PointData>\n"
             "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << off << "\"/>\n";
        off += blocksize(vel.size()*sizeof(float));
        for (int c = 0; c < ncomp; c++)
        {
            f << "        <DataArray type=\"Float32\" Name=\"" << varnames[c]
              << "\" format=\"appended\" offset=\"" << off << "\"/>\n";
            off += blocksize(vals[c].size()*sizeof(float));
        }
        f << "      </PointData>\n"
             "    </Piece>\n"
             "  </UnstructuredGrid>\n"
             "  <AppendedData encoding=\"raw\">\n_";

        append_block(f, pts);
        append_block(f, conn);
        append_block(f, offs);
        append_block(f, types);
        append_block(f, vel);
        for (int c = 0; c < ncomp; c++)
            append_block(f, vals[c]);

        f << "\n  </AppendedData>\n"
             "</VTKFile>\n";
        f.close();
    }

    // ---- rank 0: .vtm index (blocks = levels) and .pvd time series
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        char vtm_rel[200];
        snprintf(vtm_rel, sizeof(vtm_rel), "%s.vtm", stepname);

        std::ofstream vtm(base + "/" + vtm_rel, std::ios::out | std::ios::trunc);
        vtm << "<?xml version=\"1.0\"?>\n"
               "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
               "  <vtkMultiBlockDataSet>\n";
        for (int lev = 0; lev < nlevs; lev++)
        {
            vtm << "    <Block name=\"Level_" << lev << "\" index=\"" << lev << "\">\n";
            for (int r = 0; r < nprocs; r++)
                vtm << "      <DataSet index=\"" << r << "\" name=\"rank" << r
                    << "\" file=\"" << stepname << "/lev" << lev << "_r"
                    << std::setw(4) << std::setfill('0') << r << std::setfill(' ') << ".vtu\"/>\n";
            vtm << "    </Block>\n";
        }
        vtm << "  </vtkMultiBlockDataSet>\n"
               "</VTKFile>\n";
        vtm.close();

        append_pvd(base + "/plt.pvd", vtm_rel, p->simtime);
    }

    amrex::ParallelDescriptor::Barrier();
}

#endif
