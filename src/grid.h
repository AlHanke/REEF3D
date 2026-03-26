/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#ifndef GRID_H_
#define GRID_H_

#include "increment.h"

#include <vector>

class lexer;
class ghostcell;

/*!
    * @brief The grid class is responsible for defining the computational grid in REEF3D, including the grid spacing.
    * It provides methods for initializing the grid, calculating grid spacing, and managing boundary conditions.
    * The class includes members for storing nodal and cell-centered coordinates, grid spacing.
    * It also includes methods for assigning margins and initializing sigma coordinates for vertical grids.
*/
class grid : virtual public increment
{
public:
    grid() = default;
    virtual ~grid() = default;

    void assign_margin();
    void sigma_coord_ini();

    void gridspacing(lexer* p, ghostcell *pgc);

    // Non-Uniform Mesh
    std::vector<double> XN, YN, ZN; // Nodal coordinates
    std::vector<double> XP, YP, ZP; // Cell center coordinates
    std::vector<double> DXN, DYN, DZN; // Nodal grid spacing
    std::vector<double> DXP, DYP, DZP; // Cell center grid spacing
    double *ZSN,*ZSP; // Sigma coordinates (z direction)
    double DXM; // Average grid spacing (all directions)
    double DXD,DYD; // Average grid spacing in x and y direction

    double *RN,*SN,*TN; // Temporary arrays
    double DRM,DSM,DTM;
    double *DRDXN,*DSDYN,*DTDZN;
    double *DRDXP,*DSDYP,*DTDZP;

    // boundary conditions
    int *IO,*IOSL; // 0: no BC, 1: inflow, 2: outflow
    int *DFBED;

    bool i_dir,j_dir,k_dir; // existance of directions
    double x_dir,y_dir,z_dir; // existance of directions

    int **gcin, **gcout; // inflow and outflow ghost cell coordinates (i,j,k) and direction
    int gcin_count, gcout_count; // number of inflow and outflow ghost cells

    // maxcoor
    double xcoormax,xcoormin,ycoormax,ycoormin,zcoormax,zcoormin;
    double maxlength;

    const int margin = 3;

    double originx,originy,originz; // physical coordinates of the rank's origin (i=0,j=0,k=0)
    double endx,endy,endz; // physical coordinates of the rank's end point (i=knox,j=knoy,k=knoz)
    double global_xmin,global_ymin,global_zmin; // global minimum coordinates
    double global_xmax,global_ymax,global_zmax; // global maximum coordinates

    int origin_i, origin_j, origin_k; // grid indices corresponding to the rank's origin (i=0,j=0,k=0)
    int knox,knoy,knoz; // local grid size (excluding margins)
    int gknox,gknoy,gknoz; // global grid size (excluding margins)

    int imin,imax,jmin,jmax,kmin,kmax,kmaxF; // grid index limits (including margins) for array access

    double dx;
};

#endif
