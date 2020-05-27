# Hydrogen Atom Example

In this example, we will show you how to use the library to compute the
density matrix for a one dimensional hydrogen atom. This example uses no input
files. Instead, the Hamiltonian is built in the main program (which is why we've
chosen such a simple system). The Hamiltonian is then turned into a sparse
matrix, and a solver is called on it.

There are three main points to this example.
1. How to link the NTPoly libraries to your code.
2. How to build a matrix triplet list, and use it to build a distributed sparse
   matrix.
3. How to call the solver routines.

## Background

In a somewhat unconventional manner, we will do all this using the real space
basis set. Of course, it is normally neither practical nor useful to compute
the full density matrix in real space, however the equations in real space are
very simple and will allow us to demonstrate the concepts.

We need to solve the following equation:
> (-1/2 d/dr2 - 1/r)x = Ex

The kinetic energy operator is based on a 5 point stencil central difference
formula. See the following link:
http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/
central-differences/

The potential energy operator is diagonal in real space.

## Code Outline

The following steps are carried out.
1. Process the input parameters.
2. Setup the solver system.
3. Construct a linear space to solve the equations on.
4. Construct the kinetic energy operator triplet list.
5. Construct the potential energy operator triplet list.
6. Construct the distributed sparse matrix from the triplet lists.
7. Solver for the density matrix.
8. Print the density matrix to file.

## Build System

See the pre-made matrix example for details. Build with something like:

Fortran Build Instructions:
```
mpif90 main.f90 -o example \
  -I../../Build/include \
  -L../../Build/lib -lNTPoly -fopenmp -lblas

```
C++ Build Instructions:
```
mpicxx main.cc -c \
  -I../../Source/CPlusPlus -I../../Source/C

mpif90 main.o -o example \
  -L../../Build/lib -lNTPolyCPP -lNTPolyWrapper -lNTPoly -fopenmp -lstdc++ \
  -lblas -lmpi_cxx

```

(for the intel compiler, build an intermediate main.o object using the
C++ compiler, and link with the fortran compiler using the flags:
-qopenmp -cxxlib -nofor_main. When using Clang, use -lc++ instead of -lstdc++ .
-lmpicxx is only necessary for openmpi, with mpich it should be omitted.).

And then run with:
```
mpirun -np 1 ./example \
--process_rows 1 --process_columns 1 --process_slices 1 \
--threshold 1e-6 --convergence_threshold 1e-5 --grid_points 100 \
--density Density.mtx

```

Setup python environment:
```
export PYTHONPATH=../../Build/python
```

Run with python:
```
mpirun -np 1 python main.py \
--process_rows 1 --process_columns 1 --process_slices 1 \
--threshold 1e-6 --convergence_threshold 1e-5 --grid_points 100 \
--density Density.mtx

```

## Triplet List

The main difference between this example and the Premade Matrix example is the
use of triplet lists. Sparse matrices, unlike dense matrices, are very slow
to write to. Instead, we construct an intermediary data structure which will
store all the entries of the matrix, and then use it to build the matrix. The
triplet list is a list of triplets:

> [row, column, value]

Where row and column are global indices of the matrix. Each MPI process can
contribute any entry to the matrix, but each entry can only be contributed by
one process. Finally, the fill from triplet list routine is called, which
takes care of mapping the triplet list on to the 3D process grid.

## Divide Up The Matrix

In this example, we decide to divide up the linear space in one dimension. Thus
each processor is responsible for a consecutive set of grid points, and
the corresponding consecutive matrix rows. Each process fills in the those
rows for the kinetic energy and potentially energy by constructing a triplet
list. Once the triplet lists have been built, the FillFromTripletList routine
is called to construct a matrix. FillFromTripletList is a very powerful routine.
There is no need for the user to worry about data distribution. As long as
each process builds a triplet list with independent matrix elements, it will
be successful.
