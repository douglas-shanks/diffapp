#include <petscksp.h>         
#include <petscsys.h>                               
#include <omp.h>
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "q1_fem.h"

/* Routines for FEM */

void Q12D_GNi(PetscScalar Xi[], PetscScalar Ni[])
{
    PetscScalar xi = Xi[0];
    PetscScalar eta = Xi[1];
    
    /*  Compute remapped basis function*/
    
    Ni[0] = 0.25*(1.0 - xi)*(1.0 - eta);
    Ni[1] = 0.25*(1.0 - xi)*(1.0 + eta);
    Ni[2] = 0.25*(1.0 + xi)*(1.0 +eta);
    Ni[3] = 0.25*(1.0 + xi)*(1.0 - eta);
}

void Q12D_GradNi(PetscScalar Xi[], PetscScalar GradNi[][4])
{
    /* Compute the gradients of the remapped basis functions */
    PetscScalar xi = Xi[0];
    PetscScalar eta = Xi[1];
    
    /* wrt x */
    GradNi[0][0] = -0.25*(1.0-eta);
    GradNi[0][1] = -0.25*(1.0+eta);
    GradNi[0][2] =  0.25*(1.0+eta);
    GradNi[0][3] =  0.25*(1.0-eta);
    
    /* wrt y */
    GradNi[1][0] = -0.25*(1.0-xi);
    GradNi[1][1] =  0.25*(1.0-xi);
    GradNi[1][2] =  0.25*(1.0+xi);
    GradNi[1][3] = -0.25*(1.0+xi);
    
}

void Q12D_GradNx(PetscScalar GradNi[][4], PetscScalar GradNx[][4], PetscScalar elem_coords[8], PetscScalar *det_J)
{
    PetscScalar J00, J01, J10, J11, J;
    PetscScalar iJ00, iJ01, iJ10, iJ11;
    PetscInt       i;
    
    J00 = J01 = J10 = J11 = 0.0;
    
    for (i = 0; i < 4; ++i){
        PetscScalar cx = elem_coords[2*i];
        PetscScalar cy = elem_coords[2*i+1];
        
        J00 = J00 + GradNi[0][i]*cx;  /* J_xx = dx/dxi */
        J01 = J01 + GradNi[0][i]*cy;  /* J_xy = dy/dxi */
        J10 = J10 + GradNi[1][i]*cx;  /* J_yx = dx/deta */
        J11 = J11 + GradNi[1][i]*cy;  /* J_yy = dy/deta */        
    }
    
    /* det Jac */
    J = (J00*J11) - (J01*J10);
    
    /* Inverse of Jac */
    iJ00 =  J11/J;
    iJ01 = -J01/J;
    iJ10 = -J10/J;
    iJ11 =  J00/J;
    
    /* Compute isoparametric mapping */
    for(i = 0;  i < 4; i++){
        GradNx[0][i] = GradNi[0][i]*iJ00 + GradNi[1][i]*iJ01;
        GradNx[1][i] = GradNi[0][i]*iJ10 + GradNi[1][i]*iJ11;
    }
    
    if (det_J != NULL) *det_J = J;
}

void ConstructGaussQuadrature(PetscScalar gp_xi[4][2], PetscScalar gp_weight[])
{
    PetscScalar third = sqrt(1.0/3.0);
    gp_xi[0][0] = -third; gp_xi[0][1] = -third;
    gp_xi[1][0] = -third; gp_xi[1][1] =  third;
    gp_xi[2][0] =  third; gp_xi[2][1] =  third;
    gp_xi[3][0] =  third; gp_xi[3][1] = -third;
    gp_weight[0] = 1.0;
    gp_weight[1] = 1.0;
    gp_weight[2] = 1.0;
    gp_weight[3] = 1.0;
}

/* Assemble local stiffness matrix */

void ElementStiffMatrixQ1(PetscScalar Ke[4][4], PetscScalar elem_coords[])
{
    PetscInt           ngp = 4;
    PetscScalar     gp_xi[4][2];
    PetscScalar     gp_weight[4];
    PetscInt           p, i, j, k;
    PetscScalar     GNi_p[2][4], GNx_p[2][4];
    PetscScalar     J_p, tildeD[2];
    PetscScalar     B[2][4];
    
    /* define quadrature rule */
    
    ConstructGaussQuadrature(gp_xi, gp_weight);
    
    /* evaluate integral */
    
    for(p = 0; p < ngp; ++p){
        Q12D_GradNi(gp_xi[p], GNi_p);
        Q12D_GradNx(GNi_p, GNx_p, elem_coords, &J_p);
        
        for(i = 0; i < 4; ++i){
            B[0][i] = GNx_p[0][i];
            B[1][i] = GNx_p[1][i];
        }
        
        tildeD[0] = gp_weight[p]*J_p;
        tildeD[1] = gp_weight[p]*J_p;
        
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 4; ++j){
                for(k = 0; k < 2; ++k){
                    Ke[i][j] = Ke[i][j] + B[k][i]*tildeD[k]*B[k][j];
                }
            }
        }
    }
}

/* Assemble local mass matrix */

void ElementMassMatrixQ1(PetscScalar Me[4][4], PetscScalar elem_coords[])
{
    PetscInt            ngp = 4;
    PetscScalar      gp_xi[4][2];
    PetscScalar      gp_weight[4];
    PetscInt            p, i, j;
    PetscScalar      Ni_p[4];
    PetscScalar      GNi_p[2][4], GNx_p[2][4];
    PetscScalar      J_p, fac;
    PetscReal         sum=0.0;
    
    /* define quadrature rule */
    ConstructGaussQuadrature(gp_xi, gp_weight);
    
    /* evaluate integral */
    
    for(p = 0; p < ngp; ++p){
        Q12D_GNi(gp_xi[p], Ni_p);
        Q12D_GradNi(gp_xi[p], GNi_p);
        Q12D_GradNx(GNi_p, GNx_p, elem_coords, &J_p);
        fac = gp_weight[p]*J_p;

        for(i = 0; i < 4; ++i){
            for(j = 0; j < 4; ++j){
                Me[i][j] = Me[i][j] - fac*( Ni_p[i]*Ni_p[j] );
            }
        }
    }
    
    /* mass lump, simplest row sum */
    for(i = 0; i < 1; ++i){
        for(j = 0; j < 4; ++j){
            sum += Me[0][j];
        }
    }
    
    for(i = 0; i < 4; ++i){
        for(j = 0; j < 4; ++j){
            if(i==j){
                Me[i][j] = sum;
            }
            else{
                Me[i][j] = 0.0;
            }
        }
    }    
}

/* initial condition */

void initial_cond(Vec xinit, PetscReal h, PetscInt ne, PetscInt Istart, PetscInt Iend)
{
    PetscReal energy, density, x, y;
    PetscReal *coords;
    PetscInt Ii, ix, m;
    
    m = Iend - Istart;
    PetscMalloc1(2*m, &coords);
    
    for(Ii = Istart, ix = 0; Ii < Iend; ++Ii, ++ix) {
        /* coords */
        x = h*((float)(Ii%(ne+1))); y = h*((float)(Ii/(ne + 1)));
        coords[2*ix] = x; coords[2*ix +1] = y;
        
        if( ((y >= 0.0 && y <= 0.1) && ( x>=0.0 && x<= 0.1)) ){
            energy =25.0, density = 0.1;
        }
        else if( ((y >= 0.2 && y <= 0.4) && ( x>= 0.2 && x<= 0.4)) ){
            energy =0.0001, density = 100.0;
        }
        else if( (( y>=0.6 && y<= 0.8 )&&( x>= 0.2 && x<= 0.4 )) ){
            energy =0.0001, density = 100.0;
        }
        else if( (( y>=0.2 && y<= 0.4 )&&( x>= 0.6 && x<=0.8 )) ){
            energy = 0.0001, density = 100.0;
        }
        else if( (( y>=0.6 && y<= 0.8 )&&( x>= 0.6 && x<= 0.8 )) ){
            energy = 0.0001, density = 100.0;
        }
        else{
            energy = 0.1, density = 0.1;
        }
        PetscScalar v = density*energy;
        PetscInt jj = Ii;
        VecSetValues(xinit, 1, &jj, &v, INSERT_VALUES);
    }
}

void system_matrix(Mat Smat, Mat Mmat, PetscReal h, PetscInt ne, PetscInt Istart, PetscInt Iend, PetscReal dt )
{
    PetscReal  x, y;
    PetscReal *coords;
    PetscInt Ii, ix, m, i, j;
    PetscScalar  SM[4][4], MM[4][4], elem_coords[8];
    
    m = Iend - Istart;
    PetscMalloc1(2*m, &coords);
    
    /* forms the element stiffness and mass matrix for the diffusion operator and coords */
    
    for(Ii = Istart, ix = 0; Ii < Iend; ++Ii, ++ix){
        j = Ii/(ne+1); i = Ii%(ne+1);
        
        /* coords */
        x = h*((float)(Ii%(ne+1))); y = h*((float)(Ii/(ne + 1)));
        coords[2*ix] = x; coords[2*ix +1] = y;
        
        if (i < ne && j < ne){
            
            PetscInt jj, ii, idx[4];
            PetscReal density = 0.1;
            
            idx[0] = Ii; idx[1] = Ii + 1; idx[2] = Ii + (ne+1) +1; idx[3] = Ii + (ne+1);
            
            elem_coords[0] = x; elem_coords[1] = y; elem_coords[2] = x+h; elem_coords[3] = y;
            elem_coords[4] = x+h; elem_coords[5] = y+h; elem_coords[6] = x; elem_coords[7] = y+h;
            
            /* Initialise element matrices */
            
            PetscMemzero(SM, sizeof(PetscScalar)*4*4); 
            PetscMemzero(MM, sizeof(PetscScalar)*4*4);
            
            /* Form element stiffness matrix */
            
            ElementStiffMatrixQ1(SM, elem_coords);
            
            /* Form element mass matrix */
            
            ElementMassMatrixQ1(MM, elem_coords);
            
            /* Define the density */
            if( (y >= 0.2 && y <= 0.4) && ( x>= 0.2 && x<= 0.4) ){
                density = 100.0;
            }
            else if( ( y>=0.6 && y<= 0.8 )&&( x>= 0.2 && x<= 0.4 ) ){
                density = 100.0;
            }
            else if( ( y>=0.2 && y<= 0.4 )&&( x>= 0.6 && x<= 0.8 ) ){
                density = 100.0;
            }
            else if( ( y>=0.6 && y<= 0.8 )&&( x>= 0.6 && x<= 0.8 ) ){
                density = 100.0;
            }
            else{
                density = 0.1;
            }
            
            /* Additional scaling for stiffness matrix*/
            
            for(ii = 0; ii < 4; ++ii){
                for(jj = 0; jj < 4; ++jj){
                    SM[ii][jj] = -dt*(1.0 / density)*SM[ii][jj];
                    MM[ii][jj] = MM[ii][jj];
                }
            }
            MatSetValues(Smat, 4, idx, 4, idx, (const PetscScalar*)SM, ADD_VALUES); 
            MatSetValues(Mmat, 4, idx, 4, idx, (const PetscScalar*)MM, ADD_VALUES);
        }
    }
}
