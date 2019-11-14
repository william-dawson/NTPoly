# Premade Matrix Example

In this example, we will show you how to use NTPoly to compute the
density matrix for a hamiltonian and overlap matrix stored in a file. This
example assumes you've already computed a hamiltonian and overlap matrix file,
and wrote them to file in a matrix market format file. This is the simplest
example provided.

There are three main points to this example.
1. How to use link NTPoly to your code.
2. How to use the sparse matrix I/O routines.
3. How to call the solver routines.

## Background

Imagine you have a program that can already build an overlap matrix and
hamiltonian. However, you currently lack a solver to get you the density matrix.
Before integrating NTPoly with your code, you might perhaps try it out
by having the codes communicate via the file system.

First, use your code to compute the required matrices. Then, write those
matrices to file using the standard Matrix Market format. This format is
described at the
[Matrix Market Exchange](http://math.nist.gov/MatrixMarket/formats.html>)
website.

## Code Outline

The following steps are carried out.
1. Process the input parameters.
2. Construct the process grid.
3. Read the Hamiltonian and Overlap matrix from file.
4. Setup the solver system.
5. Solve for the density matrix.
6. Print the density matrix to file.

## Build System

Once you have built NTPoly, in the build folder there should be a /lib
directory. The Fortran code needs to be linked against the contained libNTPoly.a
library file. The module files are contained in the /include folder.  For
example, I have built NTPoly in a folder called "Build". So the necessary
library file is in Build/lib and the necessary modules are in Build/include.
I can build this example using gfortran with the following command:

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
--hamiltonian Hamiltonian.mtx --overlap Overlap.mtx \
--number_of_electrons 10 --threshold 1e-6 \
--converge_overlap 1e-3 --converge_density 1e-5 \
--density Density.mtx

```

Note that we're using the Fortran wrapper to link, and as a result we
have to explicitly link against the C++ standard library. For Mac users,
note that the default clang installation doesn't support openmp, so you will
need to compile everything with g++.

In the build directory, there is also a /python folder, which is used for
linking against a python program. Python requires you to set the Python path
to this directory so that it knows where to look for the python module files.

Setup python environment:
```
export PYTHONPATH=../../Build/python
```

Run with python:
```
mpirun -np 1 python main.py \
--process_rows 1 --process_columns 1 --process_slices 1 \
--hamiltonian Hamiltonian.mtx --overlap Overlap.mtx \
--number_of_electrons 10 --threshold 1e-6 \
--converge_overlap 1e-3 --converge_density 1e-5 \
--density Density.mtx

```

## Construct the process grid.

NTPoly uses a 3 dimensional process cube to distribute the work during solver
routines. Before any calls are made, you must construct this grid. In general,
the rows and columns should be close to equal. Memory use grows linearly with
the number of process slices, so the number of slices should only be
increased if the performance scaling begins to degrade. On the K computer,
we find this happens at around 8192 cores.

## Read the Hamiltonian and Overlap matrix from file.

The matrices can be read from file by simply passing in the name
of the matrix file. The density matrix and inverse square root of the overlap
matrix are built by the solver, so in Fortran you don't need to call a
constructor of them. However, in Python you do need to first construct them
with the appropriate dimensions.

## Setup the solver system.

The solver requires two important parameters: the threshold and the convergence
threshold. The threshold is used to flush small values to zero. The convergence
threshold is used to determine if a calculation is finished. For the overlap
matrix calculation, the convergence criteria is based on the norm of the
matrix. For the density matrix, it is based on the energy. One likely will want
to choose different parameters for each part of the calculation.
Additionally, we set a third parameter in the solver parameters which is used
for load balancing.

## Solve for the density matrix.

First, the inverse square root of the matrix must be computed, and then
a call to a density matrix solver routine is made (in this case, using the
TRS2 algorithm). The chemical potential is automatically computed by the
density matrix routine.

## Print the density matrix to file.

It is possible to write out a matrix either to either a text or binary format.
Writing to text is slower, but is portable across systems. Writing to binary
should probably only be used to store intermediate results.

The output file can be compared to Density-Reference.mtx to verify its
correctness.
