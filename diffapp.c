/*
    diffapp: An application to solve the 2D linear diffusion problem 
                discretised with Q1 FEM using PETSc
*/

static char help[] = "Solves 2D diffusion problem with quad finite elements. \n";
                                
#include <petscksp.h>         
#include <petscsys.h>                               
#include <omp.h>
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "q1_fem.h"


/* Main driver for diffapp */

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **args)
{
    Mat                         Amat, Smat, Mmat;
    Vec                         xx, bb, xinit;
    KSP                        ksp;
    PC                          pc;
    
    PetscErrorCode      ierr;
    PetscInt                  T, m, M, Istart, Iend, end_step;
    PetscInt                  ne = 10, da_max_its = 10000;
    PetscReal               h, residual;
    PetscReal               dt = 4.0e-2;
    float                        da_eps;
    double                    start=0.0, finish=0.0;
    PetscInt                  numit, numit1, print_matlab;
    MPI_Comm            comm;
    PetscMPIInt            npe, mype;
    PetscScalar            temp_sum;
    
    /* Initialise PETSc, MPI, variables, etc. */
    
    ierr = PetscInitialize(&argc, &args, (char*)0, help);CHKERRQ(ierr);
    comm = PETSC_COMM_WORLD;
    ierr = MPI_Comm_rank(comm, &mype); CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm, &npe); CHKERRQ(ierr);
    
    int omp_get_max_threads();
    
    /* Print output to the file diffapp.out */
    
    FILE* out_file = fopen("diffapp.out","w");
    PetscPrintf(PETSC_COMM_WORLD, "================================================ \n");
    PetscPrintf(PETSC_COMM_WORLD, "diffapp version: %.16g\n", 1.01);
    PetscPrintf(PETSC_COMM_WORLD, "MPI version \n");
    PetscPrintf(PETSC_COMM_WORLD, "MPI task count: %d\n", npe);
    PetscPrintf(PETSC_COMM_WORLD, "OMP thread count: %d\n", omp_get_max_threads());
    
    /* Print output to the terminal */

    PetscFPrintf(PETSC_COMM_WORLD,out_file,"================================================ \n");
    PetscFPrintf(PETSC_COMM_WORLD,out_file, "diffapp version: %.16g\n", 1.00);
    PetscFPrintf(PETSC_COMM_WORLD,out_file, "MPI version \n");
    PetscFPrintf(PETSC_COMM_WORLD,out_file, "MPI task count: %d\n", npe);
    PetscFPrintf(PETSC_COMM_WORLD, out_file,"OMP thread count: %d\n", omp_get_max_threads());    
    
    /* Read from input file */
    {
        FILE* in_file = fopen("diffapp.in", "r");
        ierr = fscanf(in_file, "n_elements = %d\n"
                                "initial_timestep = %lf\n"
                                "end_step = %d\n"
                                "da_max_its = %d\n"
                                "da_eps = %g\n"
                                "print_matlab = %d\n",
                                &ne, &dt, &end_step, &da_max_its, &da_eps, &print_matlab);
        fclose(in_file);
    }
    
    /* Define the mesh spacing and global number of nodes */
    
    h = (1.0)/( (float) ne );
    M = (ne + 1)*(ne + 1);
    
    /* Print these quantities to the terminal and output file */
    
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Global number of elements:   %d\n", M);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "\n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Global number of elements:   %d\n", M);    
    
    /* Start the time */
    
    start = MPI_Wtime();
    
    /* Create stiffness matrix */
    
    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, M, M,
                                    18, NULL, 6, NULL, &Smat); CHKERRQ(ierr);
                                    
    /* Create mass matrix */
    
    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, M, M,
                                    18, NULL, 6, NULL, &Mmat); CHKERRQ(ierr);                                    
                                    
    /* Create system matrix */
    
    ierr = MatCreateAIJ(comm, PETSC_DECIDE, PETSC_DECIDE, M, M,
                                    18, NULL, 6, NULL, &Amat); CHKERRQ(ierr);
    
    ierr = MatGetOwnershipRange(Amat, &Istart, &Iend); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(Smat, &Istart, &Iend); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(Mmat, &Istart, &Iend); CHKERRQ(ierr);
    
    /* Define the local node count */
    
    m = Iend - Istart;
    
    /* Generate vectors for soln and rhs */
    
    ierr = VecCreate(comm, &xx); CHKERRQ(ierr);
    ierr = VecSetSizes(xx, m, M); CHKERRQ(ierr);
    ierr = VecSetFromOptions(xx); CHKERRQ(ierr);
    ierr = VecDuplicate(xx, &bb); CHKERRQ(ierr);
    ierr = VecSet(bb, 0.0); CHKERRQ(ierr);
    ierr = VecDuplicate(xx, &xinit); CHKERRQ(ierr);
    ierr = VecSet(xinit, 0.0); CHKERRQ(ierr);
    
    /* Form the initial condition and rhs */
    
    initial_cond(xinit, h, ne, Istart, Iend);
    
    /* Form the system matrix */

    system_matrix(Smat, Mmat,  h, ne, Istart, Iend, dt);
    
    /* Assemble the matrices and vectors */
    
    ierr = MatAssemblyBegin(Smat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Smat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Mmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Mmat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(xinit); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(xinit); CHKERRQ(ierr);
    
    /* Assemble the system matrix A = S + M */
    
    ierr = MatAXPY(Smat, 1.0, Mmat, DIFFERENT_NONZERO_PATTERN);
    
    /* Setup the solver. Then finish the KSP/PC setup */
    
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, da_eps, PETSC_DEFAULT, PETSC_DEFAULT, da_max_its); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, Smat, Smat); CHKERRQ(ierr);
    
    /* Construct the rhs b = M*b */
    
    ierr = MatMult(Mmat, xinit, bb); CHKERRQ(ierr);
    
    /* Solve */
    
    ierr = VecSet(xx, 0.0); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, bb, xx); CHKERRQ(ierr);
    ierr = VecCopy(xx, xinit); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &numit); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &residual); CHKERRQ(ierr);
    
    /* Output the iteration and residual for first step */
    
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Step :                       %d\n", 1);
    PetscPrintf(PETSC_COMM_WORLD, "Iteration count :            %d\n", numit);
    PetscPrintf(PETSC_COMM_WORLD, "2-norm of the residual :     %.16g\n", residual);
    
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "\n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Step :                          %d\n", 1);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Timestep :                         %g\n", dt);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Time :                             %g\n", dt*((float)1));
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Iteration count :             %d\n", numit);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Total iteration count :       %d\n", numit);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "2-norm of the residual :     %.16g\n", residual);    

    numit1 = numit;
    
    /* Step through time */
    
    for( T = 2; T <= end_step; ++T) {
        
        /* Construct the rhs b = M*b */
        ierr = MatMult(Mmat, xinit, bb); CHKERRQ(ierr);
        
        /* solve */
        ierr = VecSet(xx,0.0); CHKERRQ(ierr);
        ierr = KSPSetUp(ksp); CHKERRQ(ierr);
        ierr = KSPSolve(ksp,bb,xx); CHKERRQ(ierr);
        ierr = VecCopy(xx, xinit); CHKERRQ(ierr);
        ierr = KSPGetIterationNumber(ksp, &numit); CHKERRQ(ierr);
        ierr = KSPGetResidualNorm(ksp, &residual); CHKERRQ(ierr);
        numit1 +=numit;
        
        /* Print output to terminal and file */
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        PetscPrintf(PETSC_COMM_WORLD, "Step :                       %d\n", T);
        PetscPrintf(PETSC_COMM_WORLD, "Iteration count :            %d\n", numit);
        PetscPrintf(PETSC_COMM_WORLD, "2-norm of the residual :     %.16g\n", residual);
        
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "\n");
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "Step :                          %d\n", T);
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "Timestep :                         %g\n", dt);
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "Time :                             %g\n", dt*((float)T));
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "Iteration count :             %d\n", numit);
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "Total iteration count :       %d\n", numit1);
        PetscFPrintf(PETSC_COMM_WORLD, out_file, "2-norm of the residual :     %.16g\n", residual); 
    }
    
    /* End timer */
    
    finish = MPI_Wtime();
    
    ierr = VecSum(xx, &temp_sum);
    temp_sum = h*h*temp_sum;
    
    /* Print out solution */
    
    if(print_matlab == 1){
        PetscViewer viewer;
        
        ierr = PetscViewerASCIIOpen(comm, "solution.m", &viewer); CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
        ierr = VecView(xx,viewer); CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
    }
    
    /* Final print out when simulation is finished */
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Calculation is complete. \n");
    PetscPrintf(PETSC_COMM_WORLD, "Diffapp has finished. \n");
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Wall clock (s):     %.16g\n", (finish - start));
    PetscPrintf(PETSC_COMM_WORLD, "================================================ \n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "\n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file,  "Calculation is complete. \n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Diffapp has finished. \n");
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Final time :                             %g\n", dt*((float)T));
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "Total iteration count :       %d\n", numit1);
    PetscFPrintf(PETSC_COMM_WORLD, out_file, "2-norm of the residual :     %.16g\n", residual);  
    PetscFPrintf(PETSC_COMM_WORLD, out_file,  "Wall clock (s):     %.16g\n", (finish - start));
    PetscFPrintf(PETSC_COMM_WORLD, out_file,  "Temperature sum (normalised):     %.16g\n", (temp_sum));
    PetscFPrintf(PETSC_COMM_WORLD, out_file,  "================================================ \n");
    
    /* Free work space */
    
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    //ierr = PCDestroy(&pc); CHKERRQ(ierr);    
    ierr = VecDestroy(&xx); CHKERRQ(ierr);
    ierr = VecDestroy(&bb); CHKERRQ(ierr);
    ierr = VecDestroy(&xinit); CHKERRQ(ierr);
    ierr = MatDestroy(&Amat); CHKERRQ(ierr);
    ierr = MatDestroy(&Smat); CHKERRQ(ierr);
    ierr = MatDestroy(&Mmat); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return 0;
}
