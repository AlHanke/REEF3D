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

#ifndef PRINT_WSFLINE_X_H_
#define PRINT_WSFLINE_X_H_

#include"increment.h"
#include<iostream>
#include<fstream>
#include<vector>

class lexer;
class fdm;
class ghostcell;
class field;
class ioflow;
class wave_theory;

using namespace std;

class print_wsfline_x : public increment
{
public:
    print_wsfline_x(lexer*,fdm*,ghostcell*);
	virtual ~print_wsfline_x();

	void wsfline(lexer*, fdm*, ghostcell*,ioflow*);


private:
    /// Tile-local j of line q's y-coordinate in the CURRENTLY installed context,
    /// or false when the line does not pass through this tile. Resolved per call,
    /// never cached -- see wsf_locate.h.
    bool locate_j(lexer*, int line, int& jj) const;

    /// Collect this rank's (x, wsf) records for every line, over all levels and
    /// tiles. Replaces the old "one slot per local i" layout, which cannot hold
    /// an AMR hierarchy: under AMReX i is tile-local, so slots collide between
    /// tiles and between levels, and the total point count is no longer knox.
    void collect(lexer*, fdm*, ghostcell*);

    /// Gather every rank's records for line q onto rank 0, then sort by x and
    /// merge duplicates. Counts differ per rank and per step, hence gatherv.
    void assemble(lexer*, ghostcell*, int line);

    void sort(double*, double*, int, int);
    void remove_multientry(lexer*, double*, double*, int&);

    // Per line: this rank's contributions, then rank 0's assembled line.
    std::vector<std::vector<double>> xloc, wsf;
    std::vector<std::vector<double>> xloc_all, wsf_all;
    std::vector<int> wsfpoints;          // assembled point count per line
    std::vector<int> recvcount, recvdispl;

    int n,q;
    ofstream wsfout;

    double xcoor;
	
	wave_theory *pwave;
};

#endif

