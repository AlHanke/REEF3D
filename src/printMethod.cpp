/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include "printMethod.h"

#include "lexer.h"

#include "print_averaging.h"
#include "turbulence.h"
#include "heat.h"
#include "multiphase.h"
#include "vorticity.h"
#include "data.h"
#include "concentration.h"
#include "sediment.h"

printMethod::printMethod(lexer *p)
{
    switch (p->P10)
    {
        case 0:
            outputFormat = new vtk3D();
            break;
        case 1: default:
            outputFormat = new vtu3D();
            break;
        case 2:
            outputFormat = new vtr3D();
            break;
        case 3:
            outputFormat = new vts3D();
            break;
    }

    if(p->mpirank==0)
    {
        outputFormat->folder("CFD");
    }
}

printMethod::~printMethod()
{
    delete outputFormat;
}

void printMethod::calcVTKOffsets(lexer *p, const int numberOfPoints, const int numberOfCells, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    // Potentail problem with the idea to combine all offsets in one function without passing the objects as arguments:
    // When different types output a different amount of data.
    n=0;

    vtkOffsets[n]=0;
    ++n;

    // // time
    // if(p->P16==1)
    // {
    //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(int)+sizeof(double);
    //     ++n;
    // }

    // velocity
    vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints*3+sizeof(int);
    ++n;

    // pmean->offset_vtk(p,vtkOffsets,n);

    // // scalars
    // {
    //     // pressure
    //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //     ++n;
    //     // k and eps
    //     pturb->offset_vtk(p,vtkOffsets,n);
    //     // eddyv
    //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //     ++n;
    //     // phi
    //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //     ++n;
    //     // T
    //     pheat->offset_vtk(p,vtkOffsets,n);
    //     // Multiphase
    //     pmp->offset_vtk(p,vtkOffsets,n);
    //     // vorticity
    //     pvort->offset_vtk(p,vtkOffsets,n);
    //     // data
    //     pdata->offset_vtk(p,vtkOffsets,n);
    //     // concentration
    //     pconc->offset_vtk(p,vtkOffsets,n);
    //     // rho
    //     if(p->P24==1 && p->F300==0)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // viscosity
    //     if(p->P71==1)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // VOF
    //     if(p->P72==1)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // Fi
    //     if(p->A10==4)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // conc
    //     if(p->P26==1)
    //     { 
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // topo
    //     if(p->P27==1)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // sediment bedlaod
    //     if(p->P76==1)
    //         psed->offset_vtk_bedload(p,vtkOffsets,n);

    //     // sediment parameters 1
    //     if(p->P77==1)
    //         psed->offset_vtk_parameter1(p,vtkOffsets,n);

    //     // sediment parameters 2
    //     if(p->P78==1)
    //         psed->offset_vtk_parameter2(p,vtkOffsets,n);

    //     // bed shear stress
    //     if(p->P79>=1)
    //         psed->offset_vtk_bedshear(p,vtkOffsets,n);

    //     // test
    //     if(p->P23==1)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // elevation
    //     vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //     ++n;
    //     // solid
    //     if(p->P25==1)
    //     { 
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // floating
    //     if(p->P28==1)
    //     {
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    //     // walldist
    //     if(p->P29==1)
    //     {   
    //         vtkOffsets[n]=vtkOffsets[n-1]+sizeof(float)*numberOfPoints+sizeof(int);
    //         ++n;
    //     }
    // }
}
