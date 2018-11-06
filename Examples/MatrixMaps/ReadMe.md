# Matrix Map Example

In this example, we will show you how to use NTPoly's matrix mapping features.
Frequently, you will want to perform some operation on all of the elements
of a matrix. For example, you might want to double all the values on the
diagonal. Or you might want to remove any element in the left half of the
matrix. The matrix mapping features help facilitate this.

## Code Outline

The following steps are carried out.
1. Process the input parameters.
2. Construct the process grid.
3. Read the matrix from file.
4. Matrix is mapped.
6. Print the mapped matrix to file.

## Build System

Once you have built NTPoly, in the build folder there should be a /lib
directory. The Fortran code needs to be linked against the contained libNTPoly.a
library file. The module files are contained in the /include folder.  For
example, I have built NTPoly in a folder called "Build". So the necessary
library file is in Build/lib and the necessary modules are in Build/include.
I can build this example using gfortran with the following command:

Fortran Build Instructions:
mpif90 main.f90 -o example \
  -I../../Build/include \
  -L../../Build/lib -lNTPoly -fopenmp -llapack -lblas

C++ Build Instructions:
mpicxx main.cc -c \
  -I../../Source/CPlusPlus -I../../Source/C

mpif90 main.o -o example \
  -L../../Build/lib -lNTPolyCPP -lNTPolyWrapper -lNTPoly -fopenmp -lstdc++ \
  -llapack -lblas

(for the intel compiler, build an intermediate main.o object using the
C++ compiler, and link with the fortran compiler using the flags:
-qopenmp -cxxlib -nofor_main. When using Clang, use -lc++ instead of -lstdc++).

And then run with:
mpirun -np 1 ./example \
--process_slices 1 --input_matrix input.mtx --output_matrix output.mtx

In the build directory, there is also a /python folder, which is used for
linking against a python program. Python requires you to set the Python path
to this directory so that it knows where to look for the python module files.

Setup python environment:
export PYTHONPATH=../../Build/python

Run with python:
mpirun -np 1 python main.py \
--process_slices 1 --input_matrix input.mtx --output_matrix output.mtx

## Mapping Procedure - Fortran

To use the map, simply define a procedure which takes the specified arguments
(see the documentation of the `MatrixMapper`). In this case, we map a
procedure which takes just the row, column, and value of a triplet. It also
takes an optional parameter which is called supplementary data. This is an
array into which you can pass any and all external data this map might depend
on. We then transform the triplet values, and return true or false based on
whether we wish to filter the value.

## Mapping Procedure - C++/Python/ETC

The interface to the maps is more general for object oriented languages
(I plan to add those approach once Fortran 2003 support is sufficient). You
can instead define a class that derives from the base class `RealOperation` or
`ComplexOperation`. These classes are functors, so you can override their
function call operator to create a suitable routine.

Each object of the derived type has a member called `data` which is a triplet.
You can access the values of an element through this triplet, and modify
them accordingly.
