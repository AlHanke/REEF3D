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

#include "printMethodSeparated.h"

#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"

#include "print_averaging.h"
#include "turbulence.h"
#include "heat.h"
#include "multiphase.h"
#include "vorticity.h"
#include "data.h"
#include "concentration.h"
#include "sediment.h"


printMethodSeparated::printMethodSeparated(lexer *p) : printMethod(p), node(p), eta(p)
{
}

void printMethodSeparated::setup(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
}

int printMethodSeparated::print(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    pgc->start4a(p,a->test,1);
    pgc->start1(p,a->u,110);
    pgc->start2(p,a->v,111);
    pgc->start3(p,a->w,112);


    pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
    pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
    pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
    pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol4(p,a->phi,14);
    pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
    //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

    a->u.ggcpol(p);
    a->v.ggcpol(p);
    a->w.ggcpol(p);
    a->press.ggcpol(p);
    a->eddyv.ggcpol(p);
    a->phi.ggcpol(p);
    a->conc.ggcpol(p);
    a->ro.ggcpol(p);
    a->visc.ggcpol(p);
    a->phi.ggcpol(p);
    a->fb.ggcpol(p);
    a->fbh4.ggcpol(p);
    //a->test.ggcpol(p);
    

    pgc->gcparacox(p,a->phi,50);
    pgc->gcparacox(p,a->phi,50);

    pgc->gcparacox(p,a->topo,150);
    pgc->gcparacox(p,a->topo,150);
    
    //pgc->start4a(p,a->topo,159);

    pgc->gcperiodicx(p,a->press,4);

    outputFormat->extent(p,pgc);
    if(p->mpirank==0)
        parallelData(p,a,pgc,pmean,pturb,pheat,pmp,pvort,pdata,pconc,psed);

    int num=0;
    if(p->P15==1)
        num = p->printcount;
    if(p->P15==2)
        num = p->count;
    outputFormat->fileName(name,"CFD",num,p->mpirank+1);

    // Open File
    ofstream result;
    result.open(name, ios::binary);
    if(result.is_open())
    {
        n=0;

        vtkOffsets[n]=0;
        ++n;

        // velocity
        vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)*3+4;
        ++n;
        
        pmean->offset_vtk(p,a,pgc,result,vtkOffsets,n);

        // scalars
        {
            // pressure
            vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
            ++n;
            // k and eps
            pturb->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // eddyv
            vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
            ++n;
            // phi
            vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
            ++n;
            // T
            pheat->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // Multiphase
            pmp->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // vorticity
            pvort->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // data
            pdata->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // concentration
            pconc->offset_vtk(p,a,pgc,result,vtkOffsets,n);
            // rho
            if(p->P24==1 && p->F300==0)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // viscosity
            if(p->P71==1)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // VOF
            if(p->P72==1)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // Fi
            if(p->A10==4)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // conc
            if(p->P26==1)
            { 
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // topo
            if(p->P27==1)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // sediment bedlaod
            if(p->P76==1)
                psed->offset_vtk_bedload(p,pgc,result,vtkOffsets,n);

            // sediment parameters 1
            if(p->P77==1)
                psed->offset_vtk_parameter1(p,pgc,result,vtkOffsets,n);

            // sediment parameters 2
            if(p->P78==1)
                psed->offset_vtk_parameter2(p,pgc,result,vtkOffsets,n);

            // bed shear stress
            if(p->P79>=1)
                psed->offset_vtk_bedshear(p,pgc,result,vtkOffsets,n);

            // test
            if(p->P23==1)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // elevation
            vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
            ++n;
            // solid
            if(p->P25==1)
            { 
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // floating
            if(p->P28==1)
            {
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
            // walldist
            if(p->P29==1)
            {   
                vtkOffsets[n]=vtkOffsets[n-1]+4*(p->pointnum)+4;
                ++n;
            }
        }
        // end scalars
        outputFormat->offset(p,vtkOffsets,n);
        //---------------------------------------------

        outputFormat->beginning(p,result);
        
        n=0;
        result<<"<PointData>"<<endl;
        result<<"\t<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
        ++n;
        
        pmean->name_vtk(p,a,pgc,result,vtkOffsets,n);

        result<<"\t<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
        ++n;

        pturb->name_vtk(p,a,pgc,result,vtkOffsets,n);

        result<<"\t<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
        ++n;
        result<<"\t<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
        ++n;

        pheat->name_vtk(p,a,pgc,result,vtkOffsets,n);
        
        pmp->name_vtk(p,a,pgc,result,vtkOffsets,n);

        pvort->name_vtk(p,a,pgc,result,vtkOffsets,n);

        pdata->name_vtk(p,a,pgc,result,vtkOffsets,n);

        pconc->name_vtk(p,a,pgc,result,vtkOffsets,n);

        if(p->P24==1 && p->F300==0)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->P71==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"viscosity\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }
        
        if(p->P72==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"VOF\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->A10==4)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->P26==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"ST_conc\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->P27==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"topo\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }
        
        if(p->P76==1)
            psed->name_vtk_bedload(p,pgc,result,vtkOffsets,n);
        
        if(p->P77==1)
            psed->name_vtk_parameter1(p,pgc,result,vtkOffsets,n);

        if(p->P78==1)
            psed->name_vtk_parameter2(p,pgc,result,vtkOffsets,n);

        if(p->P79>=1)
            psed->name_vtk_bedshear(p,pgc,result,vtkOffsets,n);

        if(p->P23==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        result<<"\t<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
        ++n;

        if(p->P25==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->P28==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"floating\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }

        if(p->P29==1)
        {
            result<<"\t<DataArray type=\"Float32\" Name=\"walldist\"  format=\"appended\" offset=\""<<vtkOffsets[n]<<"\" />"<<endl;
            ++n;
        }
        result<<"</PointData>"<<endl;
        outputFormat->ending(result,vtkOffsets,n);
        //----------------------------------------------------------------------------
        result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

        //  Velocities
        iin=3*4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol1(a->u));
            result.write((char*)&ffn, sizeof (float));

            ffn=float(p->ipol2(a->v));
            result.write((char*)&ffn, sizeof (float));

            ffn=float(p->ipol3(a->w));
            result.write((char*)&ffn, sizeof (float));
        }

        //  time average flow parameters
        pmean->print_3D(p,a,pgc,result);

        //  Pressure
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol4press(a->press)-p->pressgage);
            result.write((char*)&ffn, sizeof (float));
        }

        //  turbulence
        pturb->print_3D(p,a,pgc,result);

        //  eddyv
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol4_a(a->eddyv));
            result.write((char*)&ffn, sizeof (float));
        }

        //  phi
        node.nodefill4(p,a,pgc,a->phi,eta);
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            if(p->P18==1)
                ffn=float(p->ipol4phi(a,a->phi));
            if(p->P18==2)
                ffn = float(eta(i,j,k));
            result.write((char*)&ffn, sizeof (float));
        }

        //  T
        pheat->print_3D(p,a,pgc,result);
        
        //  Multiphase
        pmp->print_3D(p,a,pgc,result);

        //  Vorticity
        pvort->print_3D(p,a,pgc,result);

        //  Data
        pdata->print_3D(p,a,pgc,result);

        //  Concentration
        pconc->print_3D(p,a,pgc,result);

        //  density
        if(p->P24==1 && p->F300==0)
        {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->ro));
                result.write((char*)&ffn, sizeof (float));
            }
        }
        
        //  viscosity
        if(p->P71==1)
        {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->visc));
                result.write((char*)&ffn, sizeof (float));
            }
        }
        
        //  VOF
        if(p->P72==1)
        {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->vof));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  Fi
        if(p->A10==4)
        {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4press(a->Fi));
                result.write((char*)&ffn, sizeof (float));
            }

        }

        if(p->P26==1)
        {
            //  conc
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->conc));
                result.write((char*)&ffn, sizeof (float));
            }
        }
        
        if(p->P27==1)
        {
            //  topo
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->topo));
                //ffn = float(-a->bed(i,j)+p->ZN[KP1]);
                result.write((char*)&ffn, sizeof (float));
            }
        }
        
        //  sediment bedload
        if(p->P76==1)
            psed->print_3D_bedload(p,pgc,result);
        
        //  sediment parameter 1
        if(p->P77==1)
            psed->print_3D_parameter1(p,pgc,result);

        //  sediment parameter 2
        if(p->P78==1)
            psed->print_3D_parameter2(p,pgc,result);

        //  bed shear stress
        if(p->P79>=1)
            psed->print_3D_bedshear(p,pgc,result);

        //  test
        if(p->P23==1)
        {
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->test));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  elevation
        iin=4*(p->pointnum)*3;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->pos_z()+0.5*p->DZN[KP]);
            result.write((char*)&ffn, sizeof (float));
        }

        if(p->P25==1)
        {
            //  solid
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->solid));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        if(p->P28==1)
        {
            //  floating
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->fb));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        if(p->P29==1)
        {
            //  walldist
            iin=4*(p->pointnum);
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->walld));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        // -----------------------
        
        outputFormat->structureWrite(p,a,result);

        result.close();

        return 0;
    }
    else
        return 1;
}

void printMethodSeparated::parallelData(lexer* p, fdm* a, ghostcell* pgc, print_averaging *pmean, turbulence *pturb, heat *pheat, multiphase *pmp, vorticity *pvort, data *pdata, concentration *pconc, sediment *psed)
{
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	outputFormat->parallelFileName(name, "CFD", num);

	ofstream result;
	result.open(name);

	if(result.is_open())
	{
		outputFormat->beginningParallel(p,result);

		result<<"<PPointData>"<<endl;
		result<<"\t<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
		
		pmean->name_pvtk(p,a,pgc,result);

		result<<"	<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;

		pturb->name_pvtk(p,a,pgc,result);

		result<<"\t<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
		result<<"	<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;

		pheat->name_pvtk(p,a,pgc,result);
		
		pmp->name_pvtk(p,a,pgc,result);

		pvort->name_pvtk(p,a,pgc,result);

		pdata->name_pvtk(p,a,pgc,result);

		pconc->name_pvtk(p,a,pgc,result);

		if(p->P24==1 && p->F300==0)
		result<<"	<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;

		if(p->P71==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"viscosity\"/>"<<endl;
		
		if(p->P72==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"VOF\"/>"<<endl;
		
		if(p->A10==4)
		result<<"	<PDataArray type=\"Float32\" Name=\"Fi\"/>"<<endl;

		if(p->P26==1)
		{
		result<<"	<PDataArray type=\"Float32\" Name=\"ST_conc\"/>"<<endl;
		}

		if(p->P27==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"topo\"/>"<<endl;
		
		if(p->P76==1)
		psed->name_pvtk_bedload(p,pgc,result);
		
		if(p->P77==1)
		psed->name_pvtk_parameter1(p,pgc,result);

		if(p->P78==1)
		psed->name_pvtk_parameter2(p,pgc,result);

		if(p->P79>=1)
		psed->name_pvtk_bedshear(p,pgc,result);

		if(p->P23==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;

		result<<"	<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;

		if(p->P25==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;

		if(p->P28==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"floating\"/>"<<endl;

		if(p->P29==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"walldist\"/>"<<endl;

		result<<"</PPointData>"<<endl;

		outputFormat->endingParallel(result,"CFD",p->M10,num);

		result.close();
	}
	else
	cout<<"Failed to open output file."<<endl;
}