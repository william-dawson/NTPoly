# Complex Matrix Example

In this example, we will demonstrate how to use NTPoly with complex matrices.
The interface for complex matrices is largely the same, however the the triplet
list and triplet types are different. Thus the main difference in use case
will come from setting up the matrix and extracting data, but not from the
actual solvers or algebra routines.

## Example Description

NTPoly is limited to computing the functions of Hermitian matrices. One
downside of this limitation is that we cannot compute the matrix exponential of
the adjacency matrix associated with directed graphs. Indeed, computing the
exponential of a non-hermitian matrix can be very challenging.

Recently, it has been proposed (see for example [1]) that one might transform
this problem into one involving Hermitian matrices by the following
transformation. Given a node pair (A,B) and associated matrix entry Mab, if an
edge exists from A to B, but not from B to A, then we add an edge of weight
complex number i to Mab. We also subtract i from the symmetric entry
Mba. For example, the adjacency matrix:

> [0, 1]
>
> [0, 0]

Would be transformed in to:

> [0, 1 + i]
>
> [1 - i, 0]

## Build System

See the premade matrix example for details. Build with something like:

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

## Program Usage

First you need to generate the input and comparison matrix. This can be done
using the `generate.py` program I provided.

```
python generate.py 512 input.mtx exp.mtx
```

And then run with:
```
mpirun -np 1 ./example \
--process_rows 1 --process_columns 1 --process_slices 1 \
--threshold 1e-6 --input_file input.mtx --exponential_file exp-nt.mtx

```

Setup python environment:
```
export PYTHONPATH=../../Build/python
```

Run with python:
```
mpirun -np 1 python main.py \
--process_rows 1 --process_columns 1 --process_slices 1 \
--threshold 1e-6 --input_file input.mtx --exponential_file exp-nt.mtx

```

Finally you can compare the two results:
```
python compare.py exp.mtx exp-nt.mtx
```

## Program Description

First, we read in the original matrix in the `input.mtx` file. We now need
to know which edges are going in only one direction. We can do that by
computing:

> Sym = (input + input.T)

Then dividing the diagonal entries of Sym by 2, and subtracting off the
original matrix. Then, we iterate over the the entries of the guide, and add
the complex number i to a new triplet list. Notice how this triplet list is
of a different type which is especially for complex values. However, when we
finally construct matrix, it is done as normal.

> [1] Guo, Krystal, and Bojan Mohar. "Hermitian adjacency matrix of digraphs
> and mixed graphs." Journal of Graph Theory 85, no. 1 (2017): 217-248.
