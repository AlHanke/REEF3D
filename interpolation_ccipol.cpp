/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"interpolation.h"
#include"field.h"
#include"lexer.h"


double interpolation::ccipol1(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
		

    /*cout<<"i: "<<i<<" j: "<<j<<" k: "<<k<<endl;
    cout<<p->mpirank<<" xp: "<<xp<<" yp: "<<yp<<" zp: "<<zp<<"  originz: "<<p->originz<<endl;
    cout<<p->mpirank<<" XN: "<<p->XN[IP1]<<" YP: "<<p->YP[JP1]<<" YP2: "<<p->YP[JP2]<<" ZP: "<<p->ZP[KP1]<<endl;*/

    // wa
    wa = (p->XN[IP1]-xp)/p->DXP[IP];
    
    if((p->XN[IP1]-xp)/p->DXP[IP]<0.0)
    {
    wa = (p->XN[IP2]-xp)/p->DXP[IP1];
    ++i;
    }
    
    if((p->XN[IP1]-xp)/p->DXP[IP]>1.0)
    {
    wa = (p->XN[IP]-xp)/p->DXP[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }
    
	//cout<<p->mpirank<<" wa: "<<wa<<" wb: "<<wb<<" wc: "<<wc<<endl;
    
    value = lint1(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol2(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posf_j(yp);
    k = p->posc_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YN[JP1]-yp)/p->DYP[JP];
    
    if((p->YN[JP1]-yp)/p->DYP[JP]<0.0)
    {
    wb = (p->YN[JP2]-yp)/p->DYP[JP1];
    ++j;
    }
    
    if((p->YN[JP1]-yp)/p->DYP[JP]>1.0)
    {
    wb = (p->YN[JP]-yp)/p->DYP[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }

    value = lint2(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol3(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZN[KP1]-zp)/p->DZP[KP];
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]<0.0)
    {
    wc = (p->ZN[KP2]-zp)/p->DZP[KP1];
    ++k;
    }
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]>1.0)
    {
    wc = (p->ZN[KP]-zp)/p->DZP[KM1];
    --k;
    }

    value = lint3(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
//cout<<" wa: "<<wa<<" DXN[IP]: "<<p->DXN[IP]<<" DXN[IP1]: "<<p->DXN[IP1]<<" DXN[IM1]: "<<p->DXN[IM1]<<endl;
//cout<<"i: "<<i<<" j: "<<j<<" k: "<<k<<" knox: "<<p->knox<<endl;


    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }

    value =  lint4(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;
    

    return value;
}

double interpolation::ccipol4phi(fdm *a,field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }

    value =  lint4phi(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4press(fdm *a,field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
	
    /*	
    if(p->mpirank==0)
    {
    cout<<"i: "<<i<<" j: "<<j<<" k: "<<k<<endl;
    cout<<p->mpirank<<" xp: "<<xp<<" yp: "<<yp<<" zp: "<<zp<<"  originz: "<<p->originz<<endl;
    cout<<p->mpirank<<" XN: "<<p->XN[IP1]<<" YP: "<<p->YP[JP1]<<" ZP: "<<p->ZP[KP1]<<endl;
    }
    */
    
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }
    
    //if(p->mpirank==0)
    //cout<<"wa: "<<wa<<" wb: "<<wb<<" wc: "<<wc<<endl<<endl;

    value =  lint4phi(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol1_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
		
        
    // wa
    wa = (p->XN[IP1]-xp)/p->DXP[IP];
    
    if((p->XN[IP1]-xp)/p->DXP[IP]<0.0)
    wa = (p->XN[IP2]-xp)/p->DXP[IP1];
    
    if((p->XN[IP1]-xp)/p->DXP[IP]>1.0)
    wa = (p->XN[IP]-xp)/p->DXP[IM1];
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol2_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posf_j(yp);
    k = p->posc_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    
    
    // wb
    wb = (p->YN[JP1]-yp)/p->DYP[JP];
    
    if((p->YN[JP1]-yp)/p->DYP[JP]<0.0)
    wb = (p->YN[JP2]-yp)/p->DYP[JP1];
    
    if((p->YN[JP1]-yp)/p->DYP[JP]>1.0)
    wb = (p->YN[JP]-yp)/p->DYP[JM1];
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol3_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    
    
    //wc
    wc = (p->ZN[KP1]-zp)/p->DZP[KP];
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]<0.0)
    wc = (p->ZN[KP2]-zp)/p->DZP[KP1];
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]>1.0)
    wc = (p->ZN[KP]-zp)/p->DZP[KM1];

    value = lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double interpolation::ccipol4_a(field& f, double xp, double yp, double zp)
{
    ii=i;
    jj=j;
    kk=k;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posc_k(zp);
	
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZP[KP1]-zp)/p->DZN[KP];
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]<0.0)
    {
    wc = (p->ZP[KP2]-zp)/p->DZN[KP1];
    ++k;
    }
    
    if((p->ZP[KP1]-zp)/p->DZN[KP]>1.0)
    {
    wc = (p->ZP[KP]-zp)/p->DZN[KM1];
    --k;
    }


    value =  lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}


