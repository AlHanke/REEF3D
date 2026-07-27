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

#include"strain.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

strain::strain(lexer *p) : gradient(p), epsi(p->F45*p->DXM)
{
}

void strain::wallf_update(lexer *p, fdm *a, ghostcell *pgc, fieldint &wallf)
{
	int n,q;
    
    wallf.setVal(0);
    
    GC4LOOP
    if((p->gcb4[p->level][n].bc==21 || p->gcb4[p->level][n].bc==6 || p->gcb4[p->level][n].bc==7))
    {
        GCB4_TILE(n);

        i = p->gcb4[p->level][n].i;
        j = p->gcb4[p->level][n].j;
        k = p->gcb4[p->level][n].k;
        
        wallf(i,j,k)=1;
    }
    GC_TILE_RESET;
    
    QGCDF4LOOP
    {
        GCDF4_TILE(q);
        i = p->gcdf4[q].i;
        j = p->gcdf4[q].j;
        k = p->gcdf4[q].k;
        
        wallf(i,j,k)=1;
    }
    GC_TILE_RESET;
}

double strain::sij(lexer *p, fdm *a, int ii, int jj)
{
	double s=0.0;

    if(ii==1 && jj==1)
        s = 2.0*pudx(p,a);

    if((ii==1 && jj==2) || (ii==2 && jj==1))
        s = pudy(p,a) + pvdx(p,a);

    if((ii==1 && jj==3) ||( ii==3 && jj==1))
        s = pudz(p,a) + pwdx(p,a);

    if(ii==2 && jj==2)
        s = pvdy(p,a);

    if((ii==2 && jj==3) || (ii==3 && jj==2))
        s = pvdz(p,a) + pwdy(p,a);

    if(ii==3 && jj==3)
        s = 2.0*pwdz(p,a);

	return 0.5*s;
}

double strain::Sij2(lexer *p, fdm *a)
{    
    symmetricStrainRateTensor(p,a->u,a->v,a->w);
    
    double s = sqrt(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);
    
    s = s*s;

	return s;
}

double strain::strainterm(lexer *p, fdm *a)
{
    return strainterm(p,a->u,a->v,a->w);
}

double strain::strainterm(lexer *p, field &u, field &v, field &w)
{    
    symmetricStrainRateTensor(p,u,v,w);
    
    double s = sqrt(2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23);

	return s;
}

void strain::symmetricStrainRateTensor(lexer *p, field &u, field &v, field &w)
{
    if(p->j_dir==1)
    {
        s11 = pudx(p,u);
        s22 = pvdy(p,v);
        s33 = pwdz(p,w);
        s12 = (pudy(p,u) + pvdx(p,v));
        s13 = (pudz(p,u) + pwdx(p,w));
        s23 = (pvdz(p,v) + pwdy(p,w));
    }
    
    if(p->j_dir==0)
    {
        s11 = pudx(p,u);
        s22 = 0.0;
        s33 = pwdz(p,w);
        s12 = 0.0;
        s13 = (pudz(p,u) + pwdx(p,w));
        s23 = 0.0;
    }
}

double strain::magSqrSd(lexer *p, fdm *a)
{
    return magSqrSd(p,a->u,a->v,a->w);
}

double strain::magSqrSd(lexer *p, field &u, field &v, field &w)
{
    symmetricStrainRateTensor(p,u,v,w);

    skewSymmetricStrainRateTensor(p,u,v,w);

    ss11 = (s11*s11 + 0.25*s12*s12 + 0.25*s13*s13);
    ss22 = (0.25*s12*s12 + s22*s22 + 0.25*s23*s23);
    ss33 = (0.25*s13*s13 + 0.25*s23*s23 + s33*s33);
    ss12 = (0.5*s11*s12 + 0.5*s12*s22 + 0.25*s13*s23);
    ss13 = (0.5*s11*s13 + 0.25*s12*s23 + 0.5*s13*s33);
    ss23 = (0.25*s12*s13 + 0.5*s22*s23 + 0.5*s23*s33);

    rr11 = -0.25*(r12*r12 + r13*r13);
    rr22 = -0.25*(r12*r12 + r23*r23);
    rr33 = -0.25*(r13*r13 + r23*r23);
    rr12 = -0.25*r13*r23;
    rr13 = 0.25*r12*r23;
    rr23 = -0.25*r12*r13;

    double IV_SR = ss11*rr11 + 2.0*ss12*rr12 + 2.0*ss13*rr13 + ss22*rr22 + 2.0*ss23*rr23 + ss33*rr33;

    double Strain = strainterm(p,u,v,w);
    double Omega = rotationterm(p,u,v,w);
    
    double Sd = ((1.0/24.0)*((pow(Strain, 2.0)*pow(Strain, 2.0)) + (pow(Omega, 2.0)*pow(Omega, 2.0)))) + ((2.0/12.0)*(pow(Strain, 2.0)*pow(Omega, 2.0))) + (2.0*IV_SR);

	return Sd;
}


double strain::strainplain(lexer *p, fdm *a)
{
    s11 = pudx(p,a);
    s22 = p->j_dir ? pvdy(p,a) : 0.0;
    s33 = pwdz(p,a);
    s12 = p->j_dir ? (pudy(p,a) + pvdx(p,a)) : 0.0;
    s13 = (pudz(p,a) + pwdx(p,a));
    s23 = p->j_dir ? (pvdz(p,a) + pwdy(p,a)) : 0.0;

    double s = 2.0*s11*s11 + 2.0*s22*s22 + 2.0*s33*s33 + s12*s12 + s13*s13 + s23*s23;

	return s;
}
