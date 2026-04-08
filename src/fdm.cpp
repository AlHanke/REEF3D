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
    // --- single shared MultiFab for all 48 fields ---
    m_mf(make_mf(p, 48)),
    // --- velocity vector (0-2) ---
    u(p, &m_mf,  0), v(p, &m_mf,  1), w(p, &m_mf,  2),
    // --- flux vector (3-5) ---
    F(p, &m_mf,  3), G(p, &m_mf,  4), H(p, &m_mf,  5),
    // --- external flux vector (6-8) ---
    Fext(p, &m_mf,  6), Gext(p, &m_mf,  7), Hext(p, &m_mf,  8),
    // --- 6DOF face fields (9-13) ---
    fbh1(p, &m_mf,  9), fbh2(p, &m_mf, 10), fbh3(p, &m_mf, 11),
    fbh4(p, &m_mf, 12), fbh5(p, &m_mf, 13),
    // --- cell-centred scalars (14-29) ---
    press(p, &m_mf, 14),
    Fi(p, &m_mf, 15),
    eddyv(p, &m_mf, 16),
    L(p, &m_mf, 17),
    ro(p, &m_mf, 18), dro(p, &m_mf, 19), visc(p, &m_mf, 20),
    phi(p, &m_mf, 21),
    conc(p, &m_mf, 22),
    test(p, &m_mf, 23),
    topo(p, &m_mf, 24), solid(p, &m_mf, 25),
    fb(p, &m_mf, 26),
    porosity(p, &m_mf, 27), porpart(p, &m_mf, 28),
    walld(p, &m_mf, 29),
    // --- integer / slice members (owning, unchanged) ---
    nodeval(p), nodeval2D(p),
    // --- PLIC fields (30-47) ---
    nX(p, &m_mf, 30), nY(p, &m_mf, 31), nZ(p, &m_mf, 32), Alpha(p, &m_mf, 33),
    phasemarker(p, &m_mf, 34),
    vof(p, &m_mf, 35),
    vof_nt(p, &m_mf, 36), vof_nb(p, &m_mf, 37),
    vof_st(p, &m_mf, 38), vof_sb(p, &m_mf, 39),
    vof_nte(p, &m_mf, 40), vof_ntw(p, &m_mf, 41),
    vof_nbe(p, &m_mf, 42), vof_nbw(p, &m_mf, 43),
    vof_ste(p, &m_mf, 44), vof_stw(p, &m_mf, 45),
    vof_sbe(p, &m_mf, 46), vof_sbw(p, &m_mf, 47),
    // --- PTF slices (declaration order) ---
    eta(p), eta_n(p), depth(p), WL(p),
    Fifsf(p), K(p),
    etaloc(p),
    P(p), Q(p),
    bed(p),
    rhsvec(p), M(p)
#else
    u(p), F(p), Fext(p),
    v(p), G(p), Gext(p),
    w(p), H(p), Hext(p),
    press(p),
    Fi(p),
    eddyv(p),
    L(p),
    ro(p), dro(p), visc(p),
    phi(p),
    vof(p), vof_nt(p), vof_nb(p), vof_st(p), vof_sb(p), phasemarker(p),
    vof_nte(p), vof_ntw(p), vof_nbe(p), vof_nbw(p),
    vof_ste(p), vof_stw(p), vof_sbe(p), vof_sbw(p),
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
    maxF = 0.0;
    maxG = 0.0;
    maxH = 0.0;

    gi = p->W20;
    gj = p->W21;
    gk = p->W22;
}