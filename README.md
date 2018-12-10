# diffapp
Using PETSc to solve the heat equation discretised with Q1 finite elements.

Compile: make diffapp (change PETSC_DIR to local PETSc installation)

Run: OMP_NUM_THREADS=1 ./diffapp -ksp_type cg -pc_type bjacobi (for example)

Cleanup: make allclean

Visualisation: A file solution.m is produced if print_matlab=1. 

Using Ocatve (for example)

```matlab
>> solution
>> u = Vec_0x55c95d13fb60_0 ;
>> surf(reshape(u,250,250)); shading interp; view(2); colorbar
```

![output](https://user-images.githubusercontent.com/15614951/49749591-9afc6800-fca0-11e8-85d2-680c86cc3e40.jpg)



