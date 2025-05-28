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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include<utility>

/**
 * @brief Initialize the level set function phi for the free surface
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * @param pgc Ghost cell object for boundary conditions
 * 
 * This function sets up the initial distribution of the level set function phi,
 * which represents the free surface in the simulation. Positive values indicate
 * fluid presence, negative values indicate air/void regions.
 */
void initialize::iniphi(lexer* p, fdm* a, ghostcell* pgc)
{
    double r;
    double phidiff, xdiff;
    p->phimean=p->F56;

    // Initialize all phi values to -1.0 (air/void phase)
    LOOP
        a->phi(i,j,k)=-1.0;
    
    // Apply ghost cell boundary conditions
    pgc->start4(p,a->phi,50);

    // F50: Rectangular box initialization
    if(p->F50_flag==1)
        LOOP
            if(p->XN[IP]>=p->F51 && p->XN[IP]<p->F54
            && p->YN[JP]>=p->F52 && p->YN[JP]<p->F55
            && p->ZN[KP]>=p->F53 && p->ZN[KP]<p->F56)
                a->phi(i,j,k)=1.0;

    // F57: Plane initialization using equation ax + by + cz < d
    if(p->F57_1>0||p->F57_2>0||p->F57_3>0||p->F57_4>0)
    {
        LOOP
            if(p->F57_1*p->XP[IP]+ p->F57_2*p->YP[JP]+ p->F57_3*p->ZP[KP] < p->F57_4)
                a->phi(i,j,k)=1.0;
    }

    // F58: Spherical fluid region initialization
    if(p->F58_4>0.0)
    {
        LOOP
        {
            // Calculate distance from sphere center
            r = sqrt( pow(p->XP[IP]-p->F58_1,2.0)+pow(p->YP[JP]-p->F58_2,2.0)+pow(p->ZP[KP]-p->F58_3,2.0));

            // If within sphere radius, set as fluid
            if(r<=p->F58_4)
                a->phi(i,j,k)=1.0;
        }
    }

    // F59: Cylindrical fluid region initialization
    if(p->F59_r>0.0)
    {
        LOOP
        {
            // Calculate radial distance from cylinder axis (in x-y plane)
            r = sqrt( pow(p->XP[IP]-p->F59_xm,2.0)+pow(p->YP[JP]-p->F59_ym,2.0));

            // If within cylinder radius and height bounds, set as fluid
            if(r<=p->F59_r && p->pos_z()>p->F59_zs && p->pos_z()<=p->F59_ze)
                a->phi(i,j,k)=1.0;
        }
    }

    // F60: Horizontal free surface initialization
    if(p->F60>-1.0e20)
    {
        LOOP
            a->phi(i,j,k)=p->F60-p->pos_z();
    }

    // F62/F63: Sloped free surface initialization
    if((p->F60>-1.0e20 || p->F56>-1.0e20) && p->F62>-1.0e-20&& p->F63>-1.0e-20 )
    {
        phidiff=p->F62-p->phimean;
        xdiff=p->xcoormax-p->F63;

        LOOP
            if(p->pos_x() > p->F63)
                a->phi(i,j,k)=(phidiff/xdiff)*(p->pos_x()-p->F63) + p->phimean - p->pos_z() ;
    }

    // Set mean phi value
    p->phimean=p->F56;

    // Update boundary values based on initialization
    if(p->F60>-1.0e20)
    {
        p->phimean=p->F60;
        p->phiout=p->F60;
        p->fsfin=p->F60;
        p->fsfout=p->F60;
    }
    
    if(p->F61>-1.0e20)
    {
        p->phiin=p->F61;
        p->fsfin=p->F61;
    }

    if(p->F62>-1.0e20)
    {
        p->phiout=p->F62;
        p->fsfout=p->F62;
    }
}

/**
 * @brief Initialize phi values at inlet/outlet boundaries
 * @param a FDM object containing field variables
 * @param p Lexer object containing simulation parameters
 * @param pgc Ghost cell object for boundary conditions
 * 
 * Sets phi values at inlet boundaries to maintain consistent free surface levels
 */
void initialize::iniphi_io(fdm*a, lexer* p, ghostcell* pgc)
{
    // F61: Inlet boundary phi initialization
    if(p->F61>-1.0e20)
        GC4LOOP
            if(p->gcb4[n][4]==1)  // Inlet boundary (type 1)
            {
                i=p->gcb4[n][0];
                j=p->gcb4[n][1];
                k=p->gcb4[n][2];

                // Set phi values in ghost cells based on inlet level
                a->phi(i-1,j,k)=p->F61-p->pos_z();
                a->phi(i-2,j,k)=p->F61-p->pos_z();
                a->phi(i-3,j,k)=p->F61-p->pos_z();
            }

    // F62: Outlet boundary phi initialization (commented out)
    /*
    if(p->F62>-1.0e20)
        GC4LOOP
            if(p->gcb4[n][4]==2)  // Outlet boundary (type 2)
            {
                i=p->gcb4[n][0];
                j=p->gcb4[n][1];
                k=p->gcb4[n][2];

                a->phi(i+1,j,k)=p->F62-p->pos_z();
                a->phi(i+2,j,k)=p->F62-p->pos_z();
                a->phi(i+3,j,k)=p->F62-p->pos_z();
            }
    */
}

/**
 * @brief Initialize phi using box-shaped regions
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * @param pgc Ghost cell object for boundary conditions
 * 
 * Allows definition of multiple rectangular regions with different phi values
 */
void initialize::iniphi_box(lexer* p, fdm *a, ghostcell* pgc)
{
    int istart, iend, jstart, jend, kstart, kend;
    int qn;
    
    // Initialize all to air if F70 boxes are defined
    if(p->F70>0)
        LOOP
            a->phi(i,j,k)=-1.0;

    // F70: Fluid boxes (phi = 1)
    for(qn=0;qn<p->F70;++qn)
    {
        // Convert coordinates to grid indices
        istart = p->posc_i(p->F70_xs[qn]);
        iend = p->posc_i(p->F70_xe[qn]);
        
        jstart = p->posc_j(p->F70_ys[qn]);
        jend = p->posc_j(p->F70_ye[qn]);
        
        kstart = p->posc_k(p->F70_zs[qn]);
        kend = p->posc_k(p->F70_ze[qn]);

        // Set fluid phase in box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                a->phi(i,j,k)=1;
    }
    
    // F71: Air boxes (phi = -1)
    for(qn=0;qn<p->F71;++qn)
    {
        // Convert coordinates to grid indices
        istart = p->posc_i(p->F71_xs[qn]);
        iend = p->posc_i(p->F71_xe[qn]);
        
        jstart = p->posc_j(p->F71_ys[qn]);
        jend = p->posc_j(p->F71_ye[qn]);
        
        kstart = p->posc_k(p->F71_zs[qn]);
        kend = p->posc_k(p->F71_ze[qn]);

        // Set air phase in box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                a->phi(i,j,k)=-1;
    }
    
    // F72: Horizontal surface boxes
    for(qn=0;qn<p->F72;++qn)
    {
        // Convert coordinates to grid indices (only x-y, z handled separately)
        istart = p->posc_i(p->F72_xs[qn]);
        iend = p->posc_i(p->F72_xe[qn]);
        
        jstart = p->posc_j(p->F72_ys[qn]);
        jend = p->posc_j(p->F72_ye[qn]);

        // Set level set function based on height
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend)
                a->phi(i,j,k)=p->F72_h[qn]-p->pos_z();
    }
}

/**
 * @brief Initialize phi using wedge-shaped geometries
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * 
 * Creates inclined surfaces and wedge shapes for complex initial conditions
 */
void initialize::iniphi_wedge(lexer* p, fdm* a)
{
    // F112: Wedge sloped in x-direction
    for(int qn=0; qn<p->F112; ++qn)
    {
        // Calculate slope in x-z plane
        double slope=(p->F112_ze[qn]-p->F112_zs[qn])/(p->F112_xe[qn]-p->F112_xs[qn]);
        double z = p->F112_zs[qn];
        
        // Ensure consistent orientation
        if(p->F112_ze[qn]<p->F112_zs[qn])
        {
            std::swap(p->F112_ze[qn],p->F112_zs[qn]);
            z = p->F112_ze[qn];
        }

        // Apply wedge condition
        LOOP
            if(p->pos_x()>=p->F112_xs[qn] && p->pos_x()<p->F112_xe[qn]
                && p->pos_y()>=p->F112_ys[qn] && p->pos_y()<p->F112_ye[qn]
                && p->pos_z()>=p->F112_zs[qn] && p->pos_z()<slope*(p->pos_x()-p->F112_xs[qn])+z)
            {
                a->phi(i,j,k)=1.0;
            }
    }

    // F113: Wedge sloped in y-direction
    for(int qn=0; qn<p->F113; ++qn)
    {
        // Calculate slope in y-z plane
        double slope=(p->F113_ze[qn]-p->F113_zs[qn])/(p->F113_ye[qn]-p->F113_ys[qn]);
        double z = p->F113_zs[qn];
        
        // Ensure consistent orientation
        if(p->F113_ze[qn]<p->F113_zs[qn])
        {
            std::swap(p->F113_ze[qn],p->F113_zs[qn]);
            z = p->F113_ze[qn];
        }

        // Apply wedge condition
        LOOP
            if(p->pos_x()>=p->F113_xs[qn] && p->pos_x()<p->F113_xe[qn]
                && p->pos_y()>=p->F113_ys[qn] && p->pos_y()<p->F113_ye[qn]
                && p->pos_z()>=p->F113_zs[qn] && p->pos_z()<slope*(p->pos_y()-p->F113_ys[qn])+z)
            {
                a->phi(i,j,k)=1.0;
            }
    }

    // F114: Inverted wedge sloped in x-direction (slope from top)
    for(int qn=0; qn<p->F114; ++qn)
    {
        // Calculate slope in x-z plane
        double slope=(p->F114_ze[qn]-p->F114_zs[qn])/(p->F114_xe[qn]-p->F114_xs[qn]);
        double z = p->F114_zs[qn];
        
        // Ensure consistent orientation
        if(p->F114_ze[qn]<p->F114_zs[qn])
        {
            std::swap(p->F114_ze[qn],p->F114_zs[qn]);
            z = p->F114_ze[qn];
        }

        // Apply inverted wedge condition (slope starts from bottom)
        LOOP
            if(p->pos_x()>=p->F114_xs[qn] && p->pos_x()<p->F114_xe[qn]
                && p->pos_y()>=p->F114_ys[qn] && p->pos_y()<p->F114_ye[qn]
                && p->pos_z()>=z+slope*(p->pos_x()-p->F114_xs[qn]) && p->pos_z()<p->F114_ze[qn])
            {
                a->phi(i,j,k)=1.0;
            }
    }

    // F115: Inverted wedge sloped in y-direction (slope from top)
    for(int qn=0; qn<p->F115; ++qn)
    {
        // Calculate slope in y-z plane
        double slope=(p->F115_ze[qn]-p->F115_zs[qn])/(p->F115_ye[qn]-p->F115_ys[qn]);
        double z = p->F115_zs[qn];
        
        // Ensure consistent orientation
        if(p->F115_ze[qn]<p->F115_zs[qn])
        {
            std::swap(p->F115_ze[qn],p->F115_zs[qn]);
            z = p->F115_ze[qn];
        }

        // Apply inverted wedge condition (slope starts from bottom)
        LOOP
            if(p->pos_x()>=p->F115_xs[qn] && p->pos_x()<p->F115_xe[qn]
                && p->pos_y()>=p->F115_ys[qn] && p->pos_y()<p->F115_ye[qn]
                && p->pos_z()>=z+slope*(p->pos_y()-p->F115_ys[qn]) && p->pos_z()<p->F115_ze[qn])
            {
                a->phi(i,j,k)=1.0;
            }
    }
}

/**
 * @brief Calculate the surface area of the free surface interface
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * @param pgc Ghost cell object for boundary conditions
 * 
 * Computes the total interface area using the level set function gradients
 * and a Dirac delta function approximation
 */
void initialize::iniphi_surfarea(lexer* p, fdm *a, ghostcell* pgc)
{
    double dx,dy,dz,dnorm,dirac;
    double area=0.0;
    double epsi = 1.6*p->DXM;  // Interface thickness parameter
    
    LOOP
    {
        // Adaptive interface thickness based on local grid size
        epsi = (1.6/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
            
        // Calculate phi gradients using central differences
        dx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
        dy = (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
        dz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
        
        // Calculate norm of gradient (interface normal magnitude)
        dnorm = sqrt(p->DXN[IP]*p->DXN[IP] + p->DYN[JP]*p->DYN[JP] + p->DZN[KP]*p->DZN[KP]);
        
        // Initialize Dirac delta function
        dirac=0.0;
        
        // Smooth Dirac delta approximation near interface (|phi| < epsi)
        if(fabs(a->phi(i,j,k))<epsi)
            dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
        
        // Accumulate surface area contribution
        area +=  pow(p->DXM,3.0) * dirac *dnorm;
    }
    
    // Sum area contributions across all processors
    area = pgc->globalsum(area);
}

/**
 * @brief Initialize density and viscosity fields based on phi distribution
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * @param pgc Ghost cell object for boundary conditions
 * 
 * Sets material properties (density, viscosity) based on the level set function
 * using a smooth Heaviside function transition
 */
void initialize::iniphi_fields(lexer* p, fdm* a, ghostcell* pgc)
{
    // Update ghost cells for phi
    pgc->start4(p,a->phi,50);
    
    double H;  // Heaviside function value
    
    BASELOOP
    {
        // Smooth Heaviside function for material property transition
        if(a->phi(i,j,k)>p->psi)
            H=1.0;  // Pure fluid phase
        else if(a->phi(i,j,k)<-p->psi)
            H=0.0;  // Pure air phase
        else
            // Smooth transition in interface region
            H=0.5*(1.0 + a->phi(i,j,k)/(p->psi) + (1.0/PI)*sin((PI*a->phi(i,j,k))/(p->psi)));

        // Set density: W1 (fluid density), W3 (air density)
        a->ro(i,j,k)= p->W1*H + p->W3*(1.0-H);
        
        // Set viscosity: W2 (fluid viscosity), W4 (air viscosity)
        a->visc(i,j,k)= p->W2*H + p->W4*(1.0-H);
    }

    // Update ghost cells for density and viscosity
    pgc->start4(p,a->ro,1);
    pgc->start4(p,a->visc,1);
}
