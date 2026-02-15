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
Authors: Hans Bihs, Alexander Hanke (@AlHanke)
--------------------------------------------------------------------*/

#ifndef ITERATORS1D_H_
#define ITERATORS1D_H_

#define ORIGIN_I (amrex::lbound(p->amr_cell_mfi->tilebox()).x + max_i*level)
#define ORIGIN_J (amrex::lbound(p->amr_cell_mfi->tilebox()).y + max_j*level)
#define ORIGIN_K (amrex::lbound(p->amr_cell_mfi->tilebox()).z + max_k*level)

#define ZEROIP (0 + marge + ORIGIN_I)
#define ZEROJP (0 + marge + ORIGIN_J)
#define ZEROKP (0 + marge + ORIGIN_K)

#define IP (i + marge + ORIGIN_I)
#define IP1 (i + 1 + marge + ORIGIN_I)
#define IP2 (i + 2 + marge + ORIGIN_I)
#define IP3 (i + 3 + marge + ORIGIN_I)
#define IP4 (i + 4 + marge + ORIGIN_I)
#define IP5 (i + 5 + marge + ORIGIN_I)
#define IM1 (i - 1 + marge + ORIGIN_I)
#define IM2 (i - 2 + marge + ORIGIN_I)
#define IM3 (i - 3 + marge + ORIGIN_I)
#define IM4 (i - 4 + marge + ORIGIN_I)
#define IM5 (i - 5 + marge + ORIGIN_I)

#define JP (j + marge + ORIGIN_J)
#define JP1 (j + 1 + marge + ORIGIN_J)
#define JP2 (j + 2 + marge + ORIGIN_J)
#define JP3 (j + 3 + marge + ORIGIN_J)
#define JP4 (j + 4 + marge + ORIGIN_J)
#define JP5 (j + 5 + marge + ORIGIN_J)
#define JM1 (j - 1 + marge + ORIGIN_J)
#define JM2 (j - 2 + marge + ORIGIN_J)
#define JM3 (j - 3 + marge + ORIGIN_J)
#define JM4 (j - 4 + marge + ORIGIN_J)
#define JM5 (j - 5 + marge + ORIGIN_J)

#define KP (k + marge + ORIGIN_K)
#define KP1 (k + 1 + marge + ORIGIN_K)
#define KP2 (k + 2 + marge + ORIGIN_K)
#define KP3 (k + 3 + marge + ORIGIN_K)
#define KP4 (k + 4 + marge + ORIGIN_K)
#define KP5 (k + 5 + marge + ORIGIN_K)
#define KM1 (k - 1 + marge + ORIGIN_K)
#define KM2 (k - 2 + marge + ORIGIN_K)
#define KM3 (k - 3 + marge + ORIGIN_K)
#define KM4 (k - 4 + marge + ORIGIN_K)
#define KM5 (k - 5 + marge + ORIGIN_K)
#define KM6 (k - 6 + marge + ORIGIN_K)

#define IIP (ii + marge + ORIGIN_I)
#define IIP1 (ii + 1 + marge + ORIGIN_I)

#define JJP (jj + marge + ORIGIN_J)
#define JJP1 (jj + 1 + marge + ORIGIN_J)

#define KKP (kk + marge + ORIGIN_K)
#define KKP1 (kk + 1 + marge + ORIGIN_K)

#endif
