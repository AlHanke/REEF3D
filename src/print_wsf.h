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

#ifndef PRINT_WSF_H_
#define PRINT_WSF_H_

#include "increment.h"
#include <fstream>
#include <iostream>
#include <vector>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

class print_wsf : public increment
{
public:
    print_wsf(lexer*,fdm*,ghostcell*,int);
    virtual ~print_wsf();

    void height_gauge(lexer*, fdm*, ghostcell*, field&);

private:
    /// Resolve a gauge to indices valid in the CURRENTLY installed context.
    /// Under AMReX that is the tile installed by TILE_LOOP, so the result is
    /// TILE-LOCAL and must be consumed before the loop advances — no index is
    /// ever stored across iterations, which is why no TileCtx is needed and why
    /// a regrid cannot invalidate anything. False when the gauge is outside.
    /// See wsf_locate.h for why this is not position::posc_i.
    bool locate(lexer*, int gauge, int& ii, int& jj) const;

    /// Keep `value` if it comes from a finer level than what this gauge already
    /// holds. The selection key is "produced a surface", not "covers the
    /// column": a fine patch can cover a gauge horizontally while the surface
    /// sits outside its z-extent, and picking on coverage would then report the
    /// sentinel instead of the coarse level's valid answer.
    void record(int gauge, double value);

    double *x, *y; // pointers so input location arrays
    int gauge_num;

    std::vector<int> lev;      // finest level that produced a value, -1 if none
    std::vector<int> levmax;   // scratch: lev reduced across ranks
    std::vector<double> wsf;
    int n;
    std::ofstream wsfout;

    static constexpr int fileFlushMaxCount = 100; // Maximum number of gauges to write before flushing the output file
    static constexpr int precision = 9; // Precision for outputting floating-point numbers
};

#endif
