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

#include"patchBC_2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC_2D::patchBC_fillobj(lexer *p, ghostcell *pgc)
{
    // fill BC options

    // discharge
    for(int qn=0;qn<p->B411;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B411_ID[qn])
        {
            patch[qq]->Q_flag=true;
            patch[qq]->Q=p->B411_Q[qn];

            patch[qq]->gcb_uflag=2;
        }
    }

    // pressure
    for(int qn=0;qn<p->B412;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B412_ID[qn])
        {
            patch[qq]->pressure_flag=true;
            patch[qq]->pressure=p->B412_pressBC[qn];

            patch[qq]->gcb_pressflag=2;
        }
    }

    // waterlevel
    for(int qn=0;qn<p->B413;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B413_ID[qn])
        {
            patch[qq]->waterlevel_flag=true;
            patch[qq]->waterlevel=p->B413_h[qn];

            patch[qq]->gcb_phiflag=2;
        }
    }

    // perpendicular velocity
    for(int qn=0;qn<p->B414;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B414_ID[qn])
        {
            patch[qq]->Uio_flag=true;
            patch[qq]->Uio=p->B414_Uio[qn];

            patch[qq]->gcb_uflag=2;
        }
    }

    // velocity components
    for(int qn=0;qn<p->B415;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B415_ID[qn])
        {
            patch[qq]->velcomp_flag=true;
            patch[qq]->U=p->B415_U[qn];
            patch[qq]->V=p->B415_V[qn];
            patch[qq]->W=p->B415_W[qn];

            patch[qq]->gcb_uflag=2;
        }
    }

    // inflow angle
    for(int qn=0;qn<p->B416;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B416_ID[qn])
        {
            patch[qq]->alpha=(PI/180.0)*p->B416_alpha[qn];
            patch[qq]->sinalpha=sin(patch[qq]->alpha);
            patch[qq]->cosalpha=cos(patch[qq]->alpha);
        }
    }

    // inflow normals
    for(int qn=0;qn<p->B417;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B417_ID[qn])
        {
            patch[qq]->velcomp_flag=true;
            patch[qq]->Nx=p->B417_Nx[qn];
            patch[qq]->Ny=p->B417_Ny[qn];
            patch[qq]->Nz=p->B417_Nz[qn];
        }
    }

    // pressure outflow
    for(int qn=0;qn<p->B418;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B418_ID[qn])
            patch[qq]->pio_flag=1;
    }

    // hydrograph discharge
    for(int qn=0;qn<p->B421;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B421_ID[qn])
        {
            patch[qq]->hydroQ_flag=true;
            patch[qq]->Q_flag=true;

            patch[qq]->gcb_uflag=2;

            // read hydrograph
            patchBC_hydrograph_Q_read(p,qq,patch[qq]->ID);
        }
    }

    // hydrograph waterlevel
    for(int qn=0;qn<p->B422;++qn)
    {
        for(qq=0;qq<obj_count;++qq)
        if(patch[qq]->ID==p->B422_ID[qn])
        {
            patch[qq]->hydroFSF_flag=true;
            patch[qq]->waterlevel_flag=true;

            patch[qq]->gcb_phiflag=2;

            // read hydrograph
            patchBC_hydrograph_FSF_read(p,qq,patch[qq]->ID);
        }
    }

    /*
    111 - 222
    110 - 221
    100 - 211
    101 - 212
    011 - 122
    010 - 121
    001 - 112
    000 - 111
    */

    // fill gcb_flags
    for(qq=0;qq<obj_count;++qq)
    {
        if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==2)
            patch[qq]->gcb_flag = 222;

        else if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==1)
            patch[qq]->gcb_flag = 221;

        else if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==1)
            patch[qq]->gcb_flag = 211;

        else if(patch[qq]->gcb_uflag==2 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==2)
            patch[qq]->gcb_flag = 212;

        else if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==2)
            patch[qq]->gcb_flag = 122;

        else if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==2 && patch[qq]->gcb_phiflag==1)
            patch[qq]->gcb_flag = 121;

        else if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==2)
            patch[qq]->gcb_flag = 112;

        else if(patch[qq]->gcb_uflag==1 && patch[qq]->gcb_pressflag==1 && patch[qq]->gcb_phiflag==1)
            patch[qq]->gcb_flag = 111;
    }

    // fill gcbs
    for(qq=0;qq<obj_count;++qq)
        patch[qq]->counter=0;

    // line
    int count=0;
    for(int qn=0;qn<p->B440;++qn)
    {
        istart = p->posc_i(p->B440_xs[qn]);
        iend = p->posc_i(p->B440_xe[qn]);

        jstart = p->posc_j(p->B440_ys[qn]);
        jend = p->posc_j(p->B440_ye[qn]);

        // 4
        for(n=0;n<p->gcbsl4_count;++n)
        {
            i=p->gcbsl4[n][0];
            j=p->gcbsl4[n][1];

            if(i>=istart && i<iend && j>=jstart && j<jend && p->gcbsl4[n][2]==p->B440_face[qn] && p->gcbsl4[n][3]==21)
            {
                for(qq=0;qq<obj_count;++qq)
                if(patch[qq]->ID==p->B440_ID[qn])
                {
                    patch[qq]->gcb[patch[qq]->counter][0]=i;
                    patch[qq]->gcb[patch[qq]->counter][1]=j;
                    patch[qq]->gcb[patch[qq]->counter][3]=p->B440_face[qn];
                    ++patch[qq]->counter;

                    // convert gcb
                    p->gcbsl4[n][3]=patch[qq]->gcb_flag;
                }
            }
        }

        // 1
        for(n=0;n<p->gcbsl1_count;++n)
        {
            i=p->gcbsl1[n][0];
            j=p->gcbsl1[n][1];

            if(i>=istart && i<iend && j>=jstart && j<jend && p->gcbsl1[n][2]==p->B440_face[qn] && p->gcbsl1[n][3]==21)
            {
                for(qq=0;qq<obj_count;++qq)
                if(patch[qq]->ID==p->B440_ID[qn])
                {
                    // convert gcb
                    p->gcbsl1[n][3]=patch[qq]->gcb_flag;
                }
            }
        }

        // 2
        for(n=0;n<p->gcbsl2_count;++n)
        {
            i=p->gcbsl2[n][0];
            j=p->gcbsl2[n][1];

            if(i>=istart && i<iend && j>=jstart && j<jend && p->gcbsl2[n][2]==p->B440_face[qn] && p->gcbsl2[n][3]==21)
            {
                for(qq=0;qq<obj_count;++qq)
                if(patch[qq]->ID==p->B440_ID[qn])
                {
                    // convert gcb
                    p->gcbsl2[n][3]=patch[qq]->gcb_flag;
                }
            }
        }
    }
}
