/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"multiphase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reini.h"

/**
 * @brief Initialize multiphase level set functions for three-phase flow simulations
 * @param p Lexer object containing simulation parameters
 * @param a FDM object containing field variables
 * @param pgc Ghost cell object for boundary conditions
 * @param pflow IO flow object for inlet/outlet handling
 * @param pprint Printer object for output
 * @param pconvec Convection object for advection schemes
 * @param psolv Solver object for linear systems
 * 
 * This function initializes two level set functions (ls1, ls2) that define
 * the interfaces between three phases in multiphase flow simulations.
 * The combination of ls1 and ls2 determines which phase each cell contains.
 */
void multiphase_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, printer *pprint, convection *pconvec, solver *psolv)
{
    int istart, iend, jstart, jend, kstart, kend;
    int qn;
    double xc,yc,zc,r;
    
    // Initialize both level set functions to -1.0 (background phase)
    LOOP
    {
        ls1(i,j,k)=-1.0;
        ls2(i,j,k)=-1.0;
    }
    
    // Apply boundary conditions for both level set functions
    pgc->start4(p,ls1,50);
    pgc->start4(p,ls2,50);
        
    //=================================================================
    // LEVEL SET 1 (LS1) INITIALIZATION
    //=================================================================
    
    // F360: Vertical interface in x-direction (plane perpendicular to x-axis)
    if(p->F360>-1.0e20)
    {
        LOOP
            ls1(i,j,k)=p->F360-p->pos_x();
    }
    
    // F361: Vertical interface in y-direction (plane perpendicular to y-axis)
    if(p->F361>-1.0e20)
    {
        LOOP
            ls1(i,j,k)=p->F361-p->pos_y();
    }
    
    // F362: Horizontal interface in z-direction (plane perpendicular to z-axis)
    if(p->F362>-1.0e20)
    {
        LOOP
            ls1(i,j,k)=p->F362-p->pos_z();
    }
    
    // F363: Wedge-shaped regions sloped in z-direction for LS1
    for(int qn=0; qn<p->F363; ++qn)
    {
        ini_wedge(p,ls1,p->F363_xs[qn],p->F363_xe[qn],p->F363_ys[qn],p->F363_ye[qn],p->F363_zs[qn],p->F363_ze[qn]);
    }

    // F364: Inv. wedge-shaped regions sloped in x-z for LS1
    for(int qn=0; qn<p->F364; ++qn)
    {
        ini_wedge_inv(p,ls1,p->F364_xs[qn],p->F364_xe[qn],p->F364_ys[qn],p->F364_ye[qn],p->F364_zs[qn],p->F364_ze[qn]);
    }
    
    // F369: Tilted box regions with initial velocity for debris flow or landslides
    for(qn=0;qn<p->F369;++qn)
    {
        double xp1,zp1,xp2,zp2,xp3,zp3,xp4,zp4,x0,z0;
        double s,ls,alpha;
        double xc,zc;
        double xr,zr;
        double vel;
        
        // Extract tilted box parameters
        x0 = p->F369_x[qn];        // Origin x-coordinate
        z0 = p->F369_z[qn];        // Origin z-coordinate
        alpha = fabs(p->F369_a[qn]*(PI/180.0));  // Tilt angle in radians
        s = p->F369_s[qn];         // Box height (perpendicular to slope)
        ls = p->F369_l[qn];        // Box length (along slope)
        vel = p->F369_v[qn];       // Initial velocity magnitude
        
        // Calculate corner points of tilted rectangular box
        xp1 = x0;
        zp1 = z0;
        
        xp2 = s * cos(alpha) + x0;
        zp2 = s * sin(alpha) + z0;
        
        xp4 = ls * cos(PI-alpha) + x0;
        zp4 = ls * sin(PI-alpha) + z0;
        
        xp3 = s * cos(alpha) + xp4;
        zp3 = s * sin(alpha) + zp4;
        
        // Check each grid point to see if it's inside the tilted box
        LOOP
        {
            xc = p->pos_x();
            zc = p->pos_z();
            
            // Use geometric line equations to determine if point is inside quadrilateral
            // g1: Line from P1 to P2
            xr = fz(xp1,zp1,xp2,zp2,zc);
            zr = fx(xp1,zp1,xp2,zp2,xc);
            
            if(xc<xr && zc>zr)
            {
                // g2: Line from P4 to P3
                xr = fz(xp4,zp4,xp3,zp3,zc);
                zr = fx(xp4,zp4,xp3,zp3,xc);
                
                if(xc>xr && zc<zr)
                {
                    // g3: Line from P3 to P2
                    xr = fz(xp3,zp3,xp2,zp2,zc);
                    zr = fx(xp3,zp3,xp2,zp2,xc);
                    
                    if(xc<xr && zc<zr)
                    {
                        // g4: Line from P4 to P1
                        xr = fz(xp4,zp4,xp1,zp1,zc);
                        zr = fx(xp4,zp4,xp1,zp1,xc);
                        
                        if(xc>xr && zc>zr)
                        {
                            // Point is inside tilted box - set as phase 1
                            ls1(i,j,k)=1;
                            
                            // Set initial velocity components based on tilt angle
                            a->u(i,j,k) = cos(alpha)*vel;   // Horizontal velocity
                            a->w(i,j,k) = -sin(alpha)*vel;  // Vertical velocity (downward)
                        }
                    }
                }
            }
        }
    }
    
    // F370: Rectangular box regions (phase 1) for LS1
    for(qn=0;qn<p->F370;++qn)
    {
        // Convert physical coordinates to grid indices
        istart = p->posc_i(p->F370_xs[qn]);
        iend = p->posc_i(p->F370_xe[qn]);
        
        jstart = p->posc_j(p->F370_ys[qn]);
        jend = p->posc_j(p->F370_ye[qn]);
        
        kstart = p->posc_k(p->F370_zs[qn]);
        kend = p->posc_k(p->F370_ze[qn]);
    
        // Set level set to 1 (phase 1) inside box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                ls1(i,j,k)=1;
    }
    
    // F371: Rectangular void regions (background phase) for LS1
    for(qn=0;qn<p->F371;++qn)
    {
        // Convert physical coordinates to grid indices
        istart = p->posc_i(p->F371_xs[qn]);
        iend = p->posc_i(p->F371_xe[qn]);
        
        jstart = p->posc_j(p->F371_ys[qn]);
        jend = p->posc_j(p->F371_ye[qn]);
        
        kstart = p->posc_k(p->F371_zs[qn]);
        kend = p->posc_k(p->F371_ze[qn]);

        // Set level set to -1 (background phase) inside box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                ls1(i,j,k)=-1.0;
    }
    
    // F374: Circular regions (phase 1) in x-z plane for LS1
    for(qn=0;qn<p->F374;++qn)
    {
        // Adjust coordinates relative to computational origin
        xc = p->F374_xc[qn] - p->originx;
        zc = p->F374_zc[qn] - p->originz;

        LOOP
        {
            // Calculate radial distance from circle center in x-z plane
            r = sqrt(pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

            // Set phase based on distance from center
            if(r<=p->F374_r[qn])
                ls1(i,j,k)=1.0;   // Inside circle - phase 1
            else
                ls1(i,j,k)=-1.0;  // Outside circle - background phase
        }
    }
    
    // F375: Circular void regions (background phase) in x-z plane for LS1
    for(qn=0;qn<p->F375;++qn)
    {
        // Adjust coordinates relative to computational origin
        xc = p->F375_xc[qn] - p->originx;
        zc = p->F375_zc[qn] - p->originz;

        LOOP
        {
            // Calculate radial distance from circle center in x-z plane
            r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

            // Set phase based on distance from center (inverted logic)
            if(r<=p->F375_r[qn])
                ls1(i,j,k)=-1.0;  // Inside circle - background phase
            else
                ls1(i,j,k)=1.0;   // Outside circle - phase 1
        }
    }
    
    //=================================================================
    // LEVEL SET 2 (LS2) INITIALIZATION
    //=================================================================
    
    // F380: Vertical interface in x-direction for LS2
    if(p->F380>-1.0e20)
    {
        LOOP
            ls2(i,j,k)=p->F380-p->pos_x();
    }
    
    // F381: Vertical interface in y-direction for LS2
    if(p->F381>-1.0e20)
    {
        LOOP
            ls2(i,j,k)=p->F381-p->pos_y();
    }
    
    // F382: Horizontal interface in z-direction for LS2
    if(p->F382>-1.0e20)
    {
        LOOP
            ls2(i,j,k)=p->F382-p->pos_z();
    }

    // F383: Wedge-shaped regions sloped in z-direction for LS2
    for(int qn=0; qn<p->F383; ++qn)
    {
        ini_wedge(p,ls2,p->F383_xs[qn],p->F383_xe[qn],p->F383_ys[qn],p->F383_ye[qn],p->F383_zs[qn],p->F383_ze[qn]);
    }

    // F390: Rectangular box regions (phase 2) for LS2
    for(qn=0;qn<p->F390;++qn)
    {
        // Convert physical coordinates to grid indices
        istart = p->posc_i(p->F390_xs[qn]);
        iend = p->posc_i(p->F390_xe[qn]);
        
        jstart = p->posc_j(p->F390_ys[qn]);
        jend = p->posc_j(p->F390_ye[qn]);
        
        kstart = p->posc_k(p->F390_zs[qn]);
        kend = p->posc_k(p->F390_ze[qn]);

        // Set level set to 1 (phase 2) inside box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                ls2(i,j,k)=1;
    }
    
    // F391: Rectangular void regions (background phase) for LS2
    for(qn=0;qn<p->F391;++qn)
    {
        // Convert physical coordinates to grid indices
        istart = p->posc_i(p->F391_xs[qn]);
        iend = p->posc_i(p->F391_xe[qn]);
        
        jstart = p->posc_j(p->F391_ys[qn]);
        jend = p->posc_j(p->F391_ye[qn]);
        
        kstart = p->posc_k(p->F391_zs[qn]);
        kend = p->posc_k(p->F391_ze[qn]);

        // Set level set to -1 (background phase) inside box region
        LOOP
            if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend)
                ls2(i,j,k)=-1;
    }
    
    // F394: Circular regions (phase 2) in x-z plane for LS2
    for(qn=0;qn<p->F394;++qn)
    {
        // Adjust coordinates relative to computational origin
        xc = p->F394_xc[qn] - p->originx;
        zc = p->F394_zc[qn] - p->originz;

        LOOP
        {
            // Calculate radial distance from circle center in x-z plane
            r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

            // Set phase based on distance from center
            if(r<=p->F394_r[qn])
                ls2(i,j,k)=1.0;   // Inside circle - phase 2
            else
                ls2(i,j,k)=-1.0;  // Outside circle - background phase
        }
    }
    
    // F395: Circular void regions (background phase) in x-z plane for LS2
    for(qn=0;qn<p->F395;++qn)
    {
        // Adjust coordinates relative to computational origin
        xc = p->F395_xc[qn] - p->originx;
        zc = p->F395_zc[qn] - p->originz;

        LOOP
        {
            // Calculate radial distance from circle center in x-z plane
            r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

            // Set phase based on distance from center (inverted logic)
            if(r<=p->F395_r[qn])
                ls2(i,j,k)=-1.0;  // Inside circle - background phase
            else
                ls2(i,j,k)=1.0;   // Outside circle - phase 2
        }
    }
    
    // F374: Additional processing - set LS2 to background phase where LS1 is positive
    // This ensures proper phase separation in multiphase regions
    for(qn=0;qn<p->F374;++qn)
    {
        // Adjust coordinates relative to computational origin
        xc = p->F374_xc[qn] - p->originx;
        zc = p->F374_zc[qn] - p->originz;

        LOOP
        {
            // Calculate radial distance from circle center in x-z plane
            r = sqrt( pow(p->XP[IP]-xc,2.0)+pow(p->ZP[KP]-zc,2.0));

            // Where phase 1 is present, ensure phase 2 is not
            if(r<=p->F374_r[qn])
                ls2(i,j,k)=-1.0;
        }
    }
    
    // Apply boundary conditions for both level set functions
    pgc->start4(p,ls1,50);
    pgc->start4(p,ls2,50);
    
    // Perform reinitialization to ensure proper signed distance function properties
    preini->start(a,p,ls1, pgc, pflow);
    PLAINLOOP
        a->phi(i,j,k) = ls1(i,j,k);
    preini->start(a,p,ls2, pgc, pflow);
    
    // Update material properties based on phase distribution
    update(p,a,pgc);
}

/**
 * @brief Calculate z-coordinate on a line given x-coordinate
 * @param x1,z1 First point coordinates
 * @param x2,z2 Second point coordinates  
 * @param x Query x-coordinate
 * @return z-coordinate on line at given x
 * 
 * Linear interpolation function for geometric calculations
 */
double multiphase_f::fx(double x1, double z1, double x2, double z2, double x)
{
    double f;
    
    // Linear equation: z = slope * (x - x1) + z1
    f = ((z2-z1)/(x2-x1))*(x-x1) + z1;
    
    return f;
}

/**
 * @brief Calculate x-coordinate on a line given z-coordinate
 * @param x1,z1 First point coordinates
 * @param x2,z2 Second point coordinates
 * @param z Query z-coordinate
 * @return x-coordinate on line at given z
 * 
 * Inverse linear interpolation function for geometric calculations
 */
double multiphase_f::fz(double x1, double z1, double x2, double z2, double z)
{
    double f;
    
    // Inverse linear equation: x = slope * (z - z1) + x1
    f = ((x2-x1)/(z2-z1))*(z-z1) + x1;
    
    return f;
}

/**
 * @brief Initialize wedge-shaped region for a level set function
 * @param p Lexer object containing simulation parameters
 * @param f Reference to level set field to modify
 * @param xs,xe x-direction bounds
 * @param ys,ye y-direction bounds  
 * @param zs,ze z-direction bounds (defines slope)
 * 
 * Creates a wedge shape sloped in x-z plane within specified y-bounds
 */
void multiphase_f::ini_wedge(lexer* p, field& f, double xs, double xe, double ys, double ye, double zs, double ze)
{
    // Calculate slope of wedge in x-z plane
    double slope=(ze-zs)/(xe-xs);
    double z = zs;
    
    // Ensure consistent orientation (swap if end is below start)
    if(ze<zs)
    {
        std::swap(ze,zs);
        z = ze;
    }

    // Apply wedge condition to all grid points
    LOOP
        if(p->pos_x()>=xs && p->pos_x()<xe      // Within x-bounds
            && p->pos_y()>=ys && p->pos_y()<ye   // Within y-bounds
            && p->pos_z()>=zs && p->pos_z()<slope*(p->pos_x()-xs)+z)  // Below sloped surface
        {
            f(i,j,k)=1.0;  // Set as active phase
        }
}

void multiphase_f::ini_wedge_inv(lexer* p, field& f, double xs, double xe, double ys, double ye, double zs, double ze)
{
    // Calculate slope of wedge in x-z plane
    double slope=(ze-zs)/(xe-xs);
    double z = zs;
    
    // Ensure consistent orientation (swap if end is below start)
    if(ze<zs)
    {
        std::swap(ze,zs);
        z = ze;
    }

    // Apply wedge condition to all grid points
    LOOP
        if(p->pos_x()>=xs && p->pos_x()<xe      // Within x-bounds
            && p->pos_y()>=ys && p->pos_y()<ye   // Within y-bounds
            && p->pos_z()>=slope*(p->pos_x()-xs)+z && p->pos_z()<ze)
        {
            f(i,j,k)=1.0;  // Set as active phase
        }
}
