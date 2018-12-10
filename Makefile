export PETSC_DIR=/home/douglas/petsc-3.10.2

CFLAGS = -O3 -march=native -funroll-loops -cpp -fopenmp
CCX = mpicc

ifdef DEBUG
CFLAGS =-O0 -g -O -Wall -Wextra -fbounds-check -fopenmp
CCX = mpicc
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

diffapp: diffapp.o q1_fem.o
	-${CCX} ${CFLAGS} -o diffapp diffapp.o q1_fem.o ${PETSC_KSP_LIB}
	${RM} diffapp.o q1_fem.o

allclean:
	${RM} *.o *.m *.out diffapp
