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
#include <string>
#include <vector>

namespace
{
    // one raw-appended binary block in a VTK XML file: uint64 byte count + payload
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

    constexpr uint8_t VTK_HIDDENCELL = 32; // vtkDataSetAttributes HIDDENCELL ghost flag

    std::string base64(const std::string& in)
    {
        static const char tbl[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        std::string out;
        out.reserve((in.size()+2)/3*4);
        for (size_t i = 0; i < in.size(); i += 3)
        {
            uint32_t v = (uint8_t)in[i] << 16;
            if (i+1 < in.size()) v |= (uint8_t)in[i+1] << 8;
            if (i+2 < in.size()) v |= (uint8_t)in[i+2];
            out += tbl[(v >> 18) & 63];
            out += tbl[(v >> 12) & 63];
            out += (i+1 < in.size()) ? tbl[(v >> 6) & 63] : '=';
            out += (i+2 < in.size()) ? tbl[v & 63] : '=';
        }
        return out;
    }
}

// Point-interpolated output with AMR level structure as a
// vtkPartitionedDataSetCollection: one rectilinear .vtr piece per AMReX box,
// grouped per level into a Partitions entry of the .vtpc index, so ParaView's
// block selection acts as level selection. Coarse cells covered by a finer
// level carry the vtkGhostType HIDDENCELL flag and are blanked automatically,
// making the all-levels view the composite field. Velocities are interpolated
// straight from the staggered faces (REEF3D convention: u(i,j,k) = east face
// of cell (i,j,k), v north, w top), scalars from the ghosted cell data,
// matching the legacy pre-AMReX point output. plane2D collapses the pseudo-2D
// y direction to a planar grid at y=0.
void printer_CFD::print_interp_amrex(lexer* p, fdm* a,
                                     const amrex::Vector<amrex::MultiFab>& plot_mfs,
                                     const amrex::Vector<std::string>& varnames,
                                     int num, bool plane2D)
{
    const int nlevs = p->nlevs;
    const int ncomp = (int)varnames.size();

    const std::string base = plane2D ? "REEF3D_CFD_VTPC2D" : "REEF3D_CFD_VTPC";
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

            // point extents of this piece (VTK structured extents are point-based)
            const int ilo = vbx.smallEnd(0), ihi = vbx.bigEnd(0)+1;
            const int klo = vbx.smallEnd(2), khi = vbx.bigEnd(2)+1;
            const int jcell = vbx.smallEnd(1);
            const int jlo = plane2D ? jcell : vbx.smallEnd(1);
            const int jhi = plane2D ? jcell : vbx.bigEnd(1)+1;
            const int npx = ihi-ilo+1, npy = jhi-jlo+1, npz = khi-klo+1;
            const int ncx = npx-1, ncy = std::max(npy-1,1), ncz = npz-1;

            // coordinate arrays; output transform assumed axis-aligned, as in the
            // 3D plotfile path which maps the prob box corners through Xout/Yout
            std::vector<double> xc(npx), yc(npy), zc(npz);
            for (int i = 0; i < npx; i++)
                xc[i] = p->coordinates::Xout(x0 + (ilo+i)*dx, y0 + jcell*dy);
            if (plane2D)
                yc[0] = 0.0;
            else
                for (int j = 0; j < npy; j++)
                    yc[j] = p->coordinates::Yout(x0 + ilo*dx, y0 + (jlo+j)*dy);
            for (int k = 0; k < npz; k++)
                zc[k] = z0 + (klo+k)*dz;

            // point data, VTK structured order: x fastest, then y, then z
            const size_t npts = (size_t)npx*npy*npz;
            std::vector<float> vel(3*npts);
            std::vector<std::vector<float>> vals(ncomp, std::vector<float>(npts));

            size_t n = 0;
            for (int k = klo; k <= khi; k++)
            for (int j = jlo; j <= jhi; j++)      // plane2D: jlo == jhi, single plane
            for (int i = ilo; i <= ihi; i++, n++)
            {
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
                vel[3*n  ] = (float)un;
                vel[3*n+1] = (float)vn;
                vel[3*n+2] = (float)wn;
                vals[0][n] = (float)un;
                vals[1][n] = (float)vn;
                vals[2][n] = (float)wn;

                for (int c = 3; c < ncomp; c++)
                {
                    double s;
                    if (c == elev_comp)
                        s = z0 + k*dz;
                    else if (plane2D)
                        s = 0.25*( cell(i-1,j,k-1,c) + cell(i,j,k-1,c)
                                 + cell(i-1,j,k  ,c) + cell(i,j,k  ,c) );
                    else
                        s = 0.125*( cell(i-1,j-1,k-1,c) + cell(i,j-1,k-1,c)
                                  + cell(i-1,j  ,k-1,c) + cell(i,j  ,k-1,c)
                                  + cell(i-1,j-1,k  ,c) + cell(i,j-1,k  ,c)
                                  + cell(i-1,j  ,k  ,c) + cell(i,j  ,k  ,c) );
                    vals[c][n] = (float)s;
                }
            }

            // cell data: blank coarse cells covered by a finer level
            std::vector<uint8_t> ghost((size_t)ncx*ncy*ncz, 0);
            if (!finest)
            {
                size_t m = 0;
                for (int k = vbx.smallEnd(2); k <= vbx.bigEnd(2); k++)
                for (int j = (plane2D ? jcell : vbx.smallEnd(1));
                     j <= (plane2D ? jcell : vbx.bigEnd(1)); j++)
                for (int i = vbx.smallEnd(0); i <= vbx.bigEnd(0); i++, m++)
                    if (mask(i,j,k) == 0)
                        ghost[m] = VTK_HIDDENCELL;
            }

            // ---- write the .vtr piece (raw-appended binary)
            char vtr_name[300];
            snprintf(vtr_name, sizeof(vtr_name), "%s/%s/lev%d_b%05d.vtr",
                     base.c_str(), stepname, lev, mfi.index());

            char extent[120];
            snprintf(extent, sizeof(extent), "%d %d %d %d %d %d",
                     ilo, ihi, jlo, jhi, klo, khi);

            std::ofstream f(vtr_name, std::ios::out | std::ios::trunc | std::ios::binary);
            f << "<?xml version=\"1.0\"?>\n"
                 "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
                 "  <RectilinearGrid WholeExtent=\"" << extent << "\">\n"
                 "    <FieldData>\n"
                 "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\"> "
              << std::scientific << p->simtime << " </DataArray>\n"
                 "    </FieldData>\n"
                 "    <Piece Extent=\"" << extent << "\">\n";

            uint64_t off = 0;
            auto blocksize = [](uint64_t nbytes){ return sizeof(uint64_t) + nbytes; };

            f << "      <PointData>\n"
                 "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << off << "\"/>\n";
            off += blocksize(vel.size()*sizeof(float));
            for (int c = 0; c < ncomp; c++)
            {
                f << "        <DataArray type=\"Float32\" Name=\"" << varnames[c]
                  << "\" format=\"appended\" offset=\"" << off << "\"/>\n";
                off += blocksize(vals[c].size()*sizeof(float));
            }
            f << "      </PointData>\n";

            f << "      <CellData>\n"
                 "        <DataArray type=\"UInt8\" Name=\"vtkGhostType\" format=\"appended\" offset=\"" << off << "\"/>\n"
                 "      </CellData>\n";
            off += blocksize(ghost.size()*sizeof(uint8_t));

            f << "      <Coordinates>\n"
                 "        <DataArray type=\"Float64\" Name=\"Xcoords\" format=\"appended\" offset=\"" << off << "\"/>\n";
            off += blocksize(xc.size()*sizeof(double));
            f << "        <DataArray type=\"Float64\" Name=\"Ycoords\" format=\"appended\" offset=\"" << off << "\"/>\n";
            off += blocksize(yc.size()*sizeof(double));
            f << "        <DataArray type=\"Float64\" Name=\"Zcoords\" format=\"appended\" offset=\"" << off << "\"/>\n";
            f << "      </Coordinates>\n"
                 "    </Piece>\n"
                 "  </RectilinearGrid>\n"
                 "  <AppendedData encoding=\"raw\">\n_";

            append_block(f, vel);
            for (int c = 0; c < ncomp; c++)
                append_block(f, vals[c]);
            append_block(f, ghost);
            append_block(f, xc);
            append_block(f, yc);
            append_block(f, zc);

            f << "\n  </AppendedData>\n"
                 "</VTKFile>\n";
            f.close();
        }
    }

    // ---- rank 0: .vtpc index (Partitions = levels, DataSets = boxes) and .pvd.
    //      The box layout (BoxArray) is replicated on all ranks, so the file list
    //      is known without communication.
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        char vtpc_rel[200];
        snprintf(vtpc_rel, sizeof(vtpc_rel), "%s.vtpc", stepname);

        std::ofstream vtpc(base + "/" + vtpc_rel, std::ios::out | std::ios::trunc);
        vtpc << "<?xml version=\"1.0\"?>\n"
                "<VTKFile type=\"vtkPartitionedDataSetCollection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
                "  <vtkPartitionedDataSetCollection>\n";
        for (int lev = 0; lev < nlevs; lev++)
        {
            vtpc << "    <Partitions index=\"" << lev << "\" name=\"Level_" << lev << "\">\n";
            const int nfabs = (int)plot_mfs[lev].boxArray().size();
            for (int b = 0; b < nfabs; b++)
            {
                char piece[200];
                snprintf(piece, sizeof(piece), "%s/lev%d_b%05d.vtr", stepname, lev, b);
                vtpc << "      <DataSet index=\"" << b << "\" file=\"" << piece << "\"/>\n";
            }
            vtpc << "    </Partitions>\n";
        }

        // the DataAssembly provides the named hierarchy ParaView uses for
        // block selection (schema per vtkXMLPartitionedDataSetCollectionWriter)
        std::string assembly = "<?xml version=\"1.0\"?>\n"
                               "<AMR type=\"vtkDataAssembly\" version=\"1.0\" id=\"0\">\n";
        for (int lev = 0; lev < nlevs; lev++)
            assembly += "  <Level_" + std::to_string(lev) + " id=\"" + std::to_string(lev+1) + "\">\n"
                        "    <dataset id=\"" + std::to_string(lev) + "\" />\n"
                        "  </Level_" + std::to_string(lev) + ">\n";
        assembly += "</AMR>\n";

        vtpc << "    <DataAssembly encoding=\"base64\">\n"
                "      " << base64(assembly) << "\n"
                "    </DataAssembly>\n";
        vtpc << "  </vtkPartitionedDataSetCollection>\n"
                "</VTKFile>\n";
        vtpc.close();

        append_pvd(base + "/plt.pvd", vtpc_rel, p->simtime);
    }

    amrex::ParallelDescriptor::Barrier();
}

#endif
