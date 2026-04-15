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
Authors: Hans Bihs, Tobias Martin, Ahmet Soydan
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"field.h"
#include"heaviside_ls.h"

void ghostcell::solid_forcing(lexer *p, fdm *a, double alpha, field& uvel, field &vvel, field& wvel,
                              field& fx, field &fy, field &fz)
{
    // Reset heaviside field
    a->fbh1.setVal(0.0, true);
    a->fbh2.setVal(0.0, true);
    a->fbh3.setVal(0.0, true);
    a->fbh4.setVal(0.0, true);

    // Calculate forcing fields
    const double uf = 0.0;
    const double vf = 0.0;
    const double wf = 0.0;

    const double inv_alpha_dt = 1.0 / (alpha * p->dt);
    const bool onlyHeaviside = p->B21==0;
    // slip & no-slip are versions of B21==1 and do not differ

    // ---------------------------------------------------
    // Construct solid heaviside function | slip & no-slip
    // ---------------------------------------------------
    FIELDLOOP_INC_EXT(
        fx,
        FIELD_MUT(fy)
        FIELD_MUT(fz)
        FIELD_MUT_MEMBER(a,fbh1)
        FIELD_MUT_MEMBER(a,fbh2)
        FIELD_MUT_MEMBER(a,fbh3)
        FIELD_MUT_MEMBER(a,fbh4)
        FIELD_MUT_MEMBER(a,fbh5)
        FIELD_CONST(uvel)
        FIELD_CONST(vvel)
        FIELD_CONST(wvel)
        FIELD_CONST_MEMBER_INC(a,solid) FIELD_CONST_MEMBER_INC(a,topo),
        FIELD_MUT_AVGDOWN(fy) FIELD_MUT_AVGDOWN(fz)
        FIELD_MUT_MEMBER_AVGDOWN(a,fbh1) FIELD_MUT_MEMBER_AVGDOWN(a,fbh2) FIELD_MUT_MEMBER_AVGDOWN(a,fbh3)
        FIELD_MUT_MEMBER_AVGDOWN(a,fbh4) FIELD_MUT_MEMBER_AVGDOWN(a,fbh5),

        // Normal vectors, heaviside, phi calculation
        double nx; double ny; double nz; double Hx; double Hy; double Hz; double H; double dirac; double phix; double phiy; double phiz;
        normal_vector_solid_topo(p, member_solid, member_topo, Hx, Hy, Hz, H, dirac, onlyHeaviside, nx, ny, nz, phix, phiy, phiz);

        // Construct the field around the solid body to adjust the tangential velocity and calculate forcing
        if(phix<=0.0)
        fx(i,j,k) += Hx*(uf - uvel(i,j,k))*inv_alpha_dt;
        else
        fx(i,j,k) += fabs(nx)*Hx*(uf - uvel(i,j,k))*inv_alpha_dt;

        if(phiy<=0.0)
        fy(i,j,k) += Hy*(vf - vvel(i,j,k))*inv_alpha_dt;
        else
        fy(i,j,k) += fabs(ny)*Hy*(vf - vvel(i,j,k))*inv_alpha_dt;

        if(phiz<=0.0)
        fz(i,j,k) += Hz*(wf - wvel(i,j,k))*inv_alpha_dt;
        else
        fz(i,j,k) += fabs(nz)*Hz*(wf - wvel(i,j,k))*inv_alpha_dt;

        member_fbh1(i,j,k) = std::min(member_fbh1(i,j,k) + Hx, 1.0);
        member_fbh2(i,j,k) = std::min(member_fbh2(i,j,k) + Hy, 1.0);
        member_fbh3(i,j,k) = std::min(member_fbh3(i,j,k) + Hz, 1.0);

        member_fbh4(i,j,k) = std::min(member_fbh4(i,j,k) + H, 1.0);

        member_fbh5(i,j,k) = 1.0-std::min(dirac,1.0);
    )

#if USE_AMREX
    // fbh1/2/3 are contiguous at comps 9-12 in a->m_mf — one batched exchange
    startBatch(p, a->m_mf, 9,
        {{static_cast<field_amrex*>(&a->fbh1), 10},
         {static_cast<field_amrex*>(&a->fbh2), 11},
         {static_cast<field_amrex*>(&a->fbh3), 12},
         {static_cast<field_amrex*>(&a->fbh4), 40}});
#else
    start1(p,a->fbh1,10);
    start2(p,a->fbh2,11);
    start3(p,a->fbh3,12);
    start4(p,a->fbh4,40);
#endif
#if USE_AMREX
    startBatch(p, *static_cast<field_amrex*>(&fx)->get_shared_mf_vec(), 0,
        {{static_cast<field_amrex*>(&fx), 10},
         {static_cast<field_amrex*>(&fy), 11},
         {static_cast<field_amrex*>(&fz), 12}});
#else
    start1(p,fx,10);
    start2(p,fy,11);
    start3(p,fz,12);
#endif
}

template<typename GenericField>
inline void ghostcell::normal_vector_solid_topo(lexer* p, GenericField& solid, GenericField& topo, double& Hx, double& Hy, double& Hz, double& H, double& dirac, bool onlyHeaviside, double& nx, double& ny, double& nz, double& phix, double& phiy, double& phiz)
{
    const double solid_ijk = solid(i,j,k);
    const double topo_ijk = topo(i,j,k);
    const double solid_ipjk = solid(i+1,j,k);
    const double solid_ijpk = solid(i,j+1,k);
    const double solid_ijkp = solid(i,j,k+1);
    const double topo_ipjk = topo(i+1,j,k);
    const double topo_ijpk = topo(i,j+1,k);
    const double topo_ijkp = topo(i,j,k+1);

    // Phi calculation – set to zero if onlyHeaviside==true
    phix = std::min(0.5*(solid_ijk + solid_ipjk), 0.5*(topo_ijk + topo_ipjk));
    phiy = std::min(0.5*(solid_ijk + solid_ijpk), 0.5*(topo_ijk + topo_ijpk));
    phiz = std::min(0.5*(solid_ijk + solid_ijkp), 0.5*(topo_ijk + topo_ijkp));

    double phival = std::min(solid_ijk, topo_ijk);

    // Heaviside function
    if(p->topoforcing==0 && p->solidread==0)
    {
        Hx = 0.0;
        Hy = 0.0;
        Hz = 0.0;

        H = 0.0;
    }
    else
    {
        double psi;
        if(!p->j_dir)
        psi = p->X41*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
        else
        psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        double phival_x;
        double phival_y;
        double phival_z;
        if(p->topoforcing>0 && p->solidread==0)
        {
            phival_x = 0.5*(topo_ijk + topo_ipjk);
            phival_y = 0.5*(topo_ijk + topo_ijpk);
            phival_z = 0.5*(topo_ijk + topo_ijkp);

            phival = topo_ijk;
        }
        else if(p->topoforcing==0 && p->solidread>0)
        {
            phival_x = 0.5*(solid_ijk + solid_ipjk);
            phival_y = 0.5*(solid_ijk + solid_ijpk);
            phival_z = 0.5*(solid_ijk + solid_ijkp);

            phival = solid_ijk;
        }
        else if(p->topoforcing>0 && p->solidread>0)
        {
            phival_x = phix;
            phival_y = phiy;
            phival_z = phiz;
        }

        Hx = heaviside_ls(-phival_x, psi);
        Hy = heaviside_ls(-phival_y, psi);
        Hz = heaviside_ls(-phival_z, psi);
        H  = heaviside_ls(-phival,   psi);
    }

    // Dirac function
    {
        double psi;
        if(!p->j_dir)
        psi = 1.1*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
        else
        psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        if(fabs(phival)<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*phival)/psi));
        else
        dirac = 0.0;
    }

    // Normal vector calculation
    if((phix <= 0.0 && phiy <= 0.0 && phiz <= 0.0) || onlyHeaviside)
    {
        nx = 0.0;
        ny = 0.0;
        nz = 0.0;

        phix = 0.0;
        phiy = 0.0;
        phiz = 0.0;
    }
    else
    {
        const double solid_imjk = solid(i-1,j,k);
        const double solid_ijmk = solid(i,j-1,k);
        const double solid_ijkm = solid(i,j,k-1);
        const double topo_imjk = topo(i-1,j,k);
        const double topo_ijmk = topo(i,j-1,k);
        const double topo_ijkm = topo(i,j,k-1);
        const double two_dx = 2*p->DXN[IP];
        const double two_dy = 2*p->DYN[JP];
        const double two_dz = 2*p->DZN[KP];

        // Solid gradient at (i,j,k) — central differences
        const double gsx = -(solid_ipjk - solid_imjk)/two_dx;
        const double gsy = -(solid_ijpk - solid_ijmk)/two_dy;
        const double gsz = -(solid_ijkp - solid_ijkm)/two_dz;
        const double norms = sqrt(gsx*gsx + gsy*gsy + gsz*gsz);
        const double inv_ns = 1.0 / (norms > 1.0e-20 ? norms : 1.0e20);

        // Topo gradient at (i,j,k) — central differences
        const double gtx = -(topo_ipjk - topo_imjk)/two_dx;
        const double gty = -(topo_ijpk - topo_ijmk)/two_dy;
        const double gtz = -(topo_ijkp - topo_ijkm)/two_dz;
        const double normt = sqrt(gtx*gtx + gty*gty + gtz*gtz);
        const double inv_nt = 1.0 / (normt > 1.0e-20 ? normt : 1.0e20);

        // Select per face: use topo gradient when solid face average is shallower
        nx = (0.5*(solid_ijk+solid_ipjk) >= 0.5*(topo_ijk+topo_ipjk)) ? gtx*inv_nt : gsx*inv_ns;
        ny = (0.5*(solid_ijk+solid_ijpk) >= 0.5*(topo_ijk+topo_ijpk)) ? gty*inv_nt : gsy*inv_ns;
        nz = (0.5*(solid_ijk+solid_ijkp) >= 0.5*(topo_ijk+topo_ijkp)) ? gtz*inv_nt : gsz*inv_ns;
    }
}
