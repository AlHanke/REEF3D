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

#ifndef REEF3D_PCH_H_
#define REEF3D_PCH_H_

// =====================================================================
// Precompiled header — third-party only.
//
// Measured with -ftime-trace over all 1334 TUs: 4986 s of compile CPU, 81 %
// of it in the frontend, and 3329 s of that just parsing headers. 2810 s of
// those headers are third-party and are re-parsed identically in ~1300 of the
// 1334 TUs. The comments below are each header's measured share.
//
// NOTHING FROM src/ BELONGS IN HERE. lexer.h, grid_amrex.h and field_amrex.h
// are edited daily; putting them in this file would make every edit a full
// 1334-TU rebuild and give the whole win back. This file should only ever
// change when a dependency is upgraded.
//
// No .cpp includes this file — the build systems force-include it
// (-include-pch in the Makefile, target_precompile_headers in CMake).
// Turn it off with `make USE_PCH=0` or -DREEF3D_USE_PCH=OFF.
// =====================================================================

#if USE_AMREX
    #include <AMReX_MFIter.H>         // 1373 s  <- src/definitions_amrex.h
    #include <AMReX_MultiFab.H>       //  519 s  <- src/looping.h
    #include <AMReX_iMultiFab.H>      //  285 s  <- src/ArrayWrapper2D.h
    #include <AMReX_FillPatchUtil.H>  //  129 s  <- src/field_amrex.h
    #include <AMReX_IntVect.H>        //   93 s  <- src/field_base.h
    #include <AMReX_MultiFabUtil.H>   //   40 s  <- src/looping.h
    #include <AMReX_BCRec.H>          //         <- src/amrex_bc_func.h
#endif

#include <Eigen/Dense>                //  157 s  <- src/6DOF_obj.h (189 TUs only;
                                  //          drop this line if the PCH gets
                                  //          unwieldy -- it is the thinnest win)

#include <fstream>                    //   76 s  <- src/sediment.h
#include <vector>                     //   64 s  <- src/net.h
#include <set>                        //   31 s  <- src/solver.h
#include <iostream>                   //   18 s  <- src/heat.h
#include <sstream>
#include <iomanip>
#include <map>
#include <algorithm>

#endif
