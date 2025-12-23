/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"hypre_struct2D.h"

#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"vec2D.h"

static void reef3d_hypre_try_enable_device_execution()
{
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_HIP) || defined(HYPRE_USING_SYCL)
    static int configured = 0;
    if(!configured)
    {
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
#if defined(HYPRE_USING_UMPIRE_UM)
#if defined(HYPRE_MEMORY_UNIFIED)
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_UNIFIED);
#endif
#else
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
#endif
        configured = 1;
    }
#endif
}

hypre_struct2D::hypre_struct2D(lexer* p,ghostcell *pgc)
{	
    reef3d_hypre_try_enable_device_execution();

    int vecsize=p->knox*p->knoy; 
    
    p->Iarray(ilower,2);
    p->Iarray(iupper,2);
    p->Darray(values,vecsize*5);
    
    make_grid(p,pgc);	  
}

hypre_struct2D::~hypre_struct2D()
{
}

void hypre_struct2D::start(lexer* p, ghostcell* pgc, slice &f, matrix2D &M, vec2D& xvec, vec2D& rhsvec, int var)
{    
	create_solvers(p,pgc);
    
    // fill for cfd
    fill_matrix(p,pgc,M,f,rhsvec);
    
    solve(p,pgc);
        
    fillbackvec(p,f,xvec,var);
	
	delete_solvers(p,pgc);
}
