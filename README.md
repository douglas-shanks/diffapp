# diffapp
We use PETSc to solve the heat equation in the unit square

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;-\nabla\left(\alpha&space;\nabla&space;u&space;\right)&space;&plus;&space;f,&space;\text{in}&space;\&space;\Omega&space;=&space;(0,1)^{2}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;-\nabla\left(\alpha&space;\nabla&space;u&space;\right)&space;&plus;&space;f,&space;\text{in}&space;\&space;\Omega&space;=&space;(0,1)^{2}." title="\frac{\partial u}{\partial t} = -\nabla\left(\alpha \nabla u \right) + f, \text{in} \ \Omega = (0,1)^{2}," /></a>

discretised with Q1 finite elements.

The purpose of this is primarily to test PETSc solver/preconditioners.

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


