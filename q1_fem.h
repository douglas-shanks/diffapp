#include <petscksp.h>         
#include <petscsys.h>                               

/* Routines for FEM */

void Q12D_GNi(PetscScalar Xi[], PetscScalar Ni[]);

void Q12D_GradNi(PetscScalar Xi[], PetscScalar GradNi[][4]);

void Q12D_GradNx(PetscScalar GradNi[][4], PetscScalar GradNx[][4], PetscScalar elem_coords[8], PetscScalar *det_J);

void ConstructGaussQuadrature(PetscScalar gp_xi[4][2], PetscScalar gp_weight[]);

void ElementStiffMatrixQ1(PetscScalar Ke[4][4], PetscScalar elem_coords[]);

void ElementMassMatrixQ1(PetscScalar Me[4][4], PetscScalar elem_coords[]);

void initial_cond(Vec xinit, PetscReal h, PetscInt ne, PetscInt Istart, PetscInt Iend);

void system_matrix(Mat Smat, Mat Mmat, PetscReal h, PetscInt ne, PetscInt Istart, PetscInt Iend, PetscReal dt );
