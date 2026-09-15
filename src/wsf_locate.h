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

#ifndef WSF_LOCATE_H_
#define WSF_LOCATE_H_

#include <vector>

class lexer;

// =====================================================================
// Shared helpers for the free-surface probes (print_wsf, print_wsfline_x,
// print_wsfline_y).
//
// The probes resolve a physical coordinate to indices ON EVERY CALL, from
// inside LEVEL_LOOP TILE_LOOP, rather than caching indices at construction.
// Three reasons, all of which a cached index gets wrong under AMReX:
//
//   * printer_CFD is built in logic_cfd(), which runs BEFORE
//     driver_ini_cfd()'s regrid loop -- a location fixed at construction only
//     ever saw a single-level grid.
//   * a tile-local (i,j,k) is meaningless without the TileCtx that produced
//     it, and nothing outlives its loop iteration here, so none is needed.
//   * an adaptive regrid moves the point between levels mid-run.
//
// Deliberately NOT built on position::posc_i/posc_j: those are broken above
// level 0 (they search a window sized by the stale level-0 `knox` and index
// XN without the level stride), and fixing them is a separate change with
// ~30 unrelated call sites.
// =====================================================================

namespace wsf_locate
{
    /// Index of the cell containing coordinate `s`, searched over the node
    /// window [org, org+imax+1] of the monotonic nodal array `N`; -1 when `s`
    /// lies outside that window.
    ///
    /// The caller passes the window as ORIGIN_* / *MAX_LOOP, which is what
    /// makes one implementation serve both builds: legacy those expand to the
    /// rank's subdomain, under AMReX to the installed tile at the installed
    /// level, with the level stride into XN/YN/ZN already folded into ORIGIN_*.
    /// The returned index is therefore relative to the same origin the field
    /// accessors use -- tile-local under AMReX.
    int cell_1d(const std::vector<double>& N, int org, int imax, double s);

    /// Whether the tile-local cell (i,j,k) at the installed level carries the
    /// authoritative value, i.e. is not overlaid by a finer level. Always true
    /// without AMReX, and on the finest level.
    ///
    /// Needed by the line probes, which emit one record per cell along the
    /// line: a coarse cell centre under a fine patch does not coincide with any
    /// fine centre (ref_ratio 2 puts them at +/- dx/4), so without this filter
    /// the coarse points survive the dedupe and interleave with the fine ones.
    /// The point probe does not need it -- it selects on the finest level that
    /// actually produced a surface, which is strictly more robust.
    bool uncovered(lexer* p, int i, int j, int k);
}

#endif
