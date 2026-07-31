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

#include"fdm.h"
#include"lexer.h"

fdm::fdm(lexer *p) :
#if USE_AMREX
    // --- single shared MultiFab for all 47 fields ---
    m_mf(make_mf(p, 47, &m_mf)),
    // --- velocity vector (0-2) ---
    u(p, &m_mf,  0), v(p, &m_mf,  1), w(p, &m_mf,  2),
    // --- flux vector (3-5) ---
    F(p, &m_mf,  3), G(p, &m_mf,  4), H(p, &m_mf,  5),
    // --- external flux vector (6-8) ---
    Fext(p, &m_mf,  6), Gext(p, &m_mf,  7), Hext(p, &m_mf,  8),
    // --- 6DOF face fields (9-13) ---
    fbh1(p, &m_mf,  9), fbh2(p, &m_mf, 10), fbh3(p, &m_mf, 11),
    fbh4(p, &m_mf, 12), fbh5(p, &m_mf, 13),
    // --- cell-centred scalars (14-28) ---
    press(p, &m_mf, 14),
    Fi(p, &m_mf, 15),
    eddyv(p, &m_mf, 16),
    L(p, &m_mf, 17),
    ro(p, &m_mf, 18), visc(p, &m_mf, 19),
    phi(p, &m_mf, 20),
    conc(p, &m_mf, 21),
    test(p, &m_mf, 22),
    topo(p, &m_mf, 23), solid(p, &m_mf, 24),
    fb(p, &m_mf, 25),
    porosity(p, &m_mf, 26), porpart(p, &m_mf, 27),
    walld(p, &m_mf, 28),
    // --- integer / slice members (owning, unchanged) ---
    nodeval(p), nodeval2D(p),
    // --- PLIC fields (29-46) ---
    nX(p, &m_mf, 29), nY(p, &m_mf, 30), nZ(p, &m_mf, 31), Alpha(p, &m_mf, 32),
    phasemarker(p, &m_mf, 33),
    vof(p, &m_mf, 34),
    vof_nt(p, &m_mf, 35), vof_nb(p, &m_mf, 36),
    vof_st(p, &m_mf, 37), vof_sb(p, &m_mf, 38),
    vof_nte(p, &m_mf, 39), vof_ntw(p, &m_mf, 40),
    vof_nbe(p, &m_mf, 41), vof_nbw(p, &m_mf, 42),
    vof_ste(p, &m_mf, 43), vof_stw(p, &m_mf, 44),
    vof_sbe(p, &m_mf, 45), vof_sbw(p, &m_mf, 46),
    // --- PTF slices (declaration order) ---
    eta(p), eta_n(p), depth(p), WL(p),
    Fifsf(p), K(p),
    etaloc(p),
    P(p), Q(p),
    bed(p),
    rhsvec(p), M(p),
    m_mf_diag(make_mf(p, 30, &m_mf_diag)),
    u0(p, &m_mf_diag, 0), du0(p, &m_mf_diag, 1), v0(p, &m_mf_diag, 2), dv0(p, &m_mf_diag, 3), w0(p, &m_mf_diag, 4), dw0(p, &m_mf_diag, 5), pcorr0(p, &m_mf_diag, 6), div0(p, &m_mf_diag, 7), phi0(p, &m_mf_diag, 8), ro0(p, &m_mf_diag, 9),
    u1(p, &m_mf_diag, 10), du1(p, &m_mf_diag, 11), v1(p, &m_mf_diag, 12), dv1(p, &m_mf_diag, 13), w1(p, &m_mf_diag, 14), dw1(p, &m_mf_diag, 15), pcorr1(p, &m_mf_diag, 16), div1(p, &m_mf_diag, 17), phi1(p, &m_mf_diag, 18), ro1(p, &m_mf_diag, 19),
    u2(p, &m_mf_diag, 20), du2(p, &m_mf_diag, 21), v2(p, &m_mf_diag, 22), dv2(p, &m_mf_diag, 23), w2(p, &m_mf_diag, 24), dw2(p, &m_mf_diag, 25), pcorr2(p, &m_mf_diag, 26), div2(p, &m_mf_diag, 27), phi2(p, &m_mf_diag, 28), ro2(p, &m_mf_diag, 29)
#else
    u(p), F(p), Fext(p),
    v(p), G(p), Gext(p),
    w(p), H(p), Hext(p),
    press(p),
    Fi(p),
    eddyv(p),
    L(p),
    ro(p), visc(p),
    phi(p),
    vof(p), vof_nt(p), vof_nb(p), vof_st(p), vof_sb(p), phasemarker(p),
    vof_nte(p), vof_ntw(p), vof_nbe(p), vof_nbw(p), vof_ste(p), vof_stw(p), vof_sbe(p), vof_sbw(p),
    conc(p),
    topo(p), solid(p),
    test(p),
    fb(p), fbh1(p), fbh2(p), fbh3(p), fbh4(p), fbh5(p),
    porosity(p), porpart(p),
    walld(p),
    nodeval(p), nodeval2D(p), etaloc(p),
    eta(p), eta_n(p), depth(p), WL(p),
    Fifsf(p), K(p),
    P(p), Q(p), bed(p),
    rhsvec(p), M(p),
    nX(p), nY(p), nZ(p), Alpha(p)
#endif
{
#if USE_AMREX
    p->register_mf(&m_mf, 47);
#endif

    maxF = 0.0;
    maxG = 0.0;
    maxH = 0.0;

    gi = p->W20;
    gj = p->W21;
    gk = p->W22;
}