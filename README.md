# diffapp
Using PETSc to solve the heat equation discretised with Q1 finite elements.

Compile:

make diffapp

Run:

OMP_NUM_THREADS=1 ./diffapp -ksp_type cg -pc_type bjacobi (for example)

