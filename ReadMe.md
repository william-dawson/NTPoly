Project Overview
================================================================================
![Build Status](https://github.com/william-dawson/NTPoly/workflows/CI/badge.svg)

NTPoly is a massively parallel library for computing the functions of sparse,
Hermitian matrices based on polynomial expansions. For sufficiently sparse
matrices, most of the matrix functions in NTPoly can be computed in linear
time.

Set Up Guide
--------------------------------------------------------------------------------
NTPoly is freely available and open source under the MIT license. It can be
downloaded from the [Github](https://github.com/william-dawson/NTPoly)
repository. We of course recommend that you download a
[release version](https://github.com/william-dawson/NTPoly/releases)
to get started.

Installing NTPoly requires the following software:

* A Fortran Compiler.
* An MPI Installation (MPI-3 Standard+).
* CMake (Version 3.2+).

The following optional software can greatly enhance the NTPoly experience:

* BLAS: for multiplying dense matrices, if they emerge in the calculation.
* A C++ Compiler for building C++ bindings.
* Ford: for building documentation.
* Doxygen: for building C++ documentation.
* SWIG (Version 3.0+): for building the Python bindings.
* Python (Version 2.7+): if you would like python bindings.
* MPI4PY: for testing.
* SciPy: for testing.
* NumPy: for testing.

NTPoly uses CMake as a build system. First, take a look in the Targets
directory. You'll find a list of `.cmake` files which have example
configurations on popular systems. You should copy one of these files, and
create your own mymachine.cmake file. Then, cd into the Build directory, and
type:
> cmake -DCMAKE_TOOLCHAIN_FILE=../Targets/mymachine.cmake ..

After that you can build using:
> make

And for the documentation:
> make doc

If you aren't cross compiling and have built the python bindings, you can
perform local tests using:
> make test

There are a few options you can pass to CMake to modify the build. A few
useful standard options are:
* `-DCMAKE_BUILD_TYPE=` `Debug` or `Release`.
* `-DCMAKE_INSTALL_PREFIX=` followed by the path to your desired install
 directory.

There are also some custom options special for NTPoly:
* `-DFORTRAN_ONLY=` set to `Yes` if you only want to build the Fortran bindings.
* `-DNOSWIG=` set to `Yes` if you don't want to build Python bindings.
* `-DUSE_MPIH=` on some systems, there is no `use mpi` feature for Fortran,
 just `#include "mpi.h"`. You can set this option to activate the later.
* `-DNOIALLGATHER=` on older MPI implementations, there is no non blocking
 collective operations. You can disable this feature using this option, but
 beware this might reduce performance.

[Online documentation](https://william-dawson.github.io/NTPoly/documentation/)
is also available. Further details about the library can be found on the
[Wiki](https://github.com/william-dawson/NTPoly/wiki).

Basic Theory
--------------------------------------------------------------------------------
The theory of matrix functions is a long studied branch of matrix algebra.
Matrix functions have a wide range of applications, including graph problems,
differential equations, and materials science. Common examples of matrix
functions include the matrix exponential:

> f(A) = e^A.

from the study of networks, or the inverse square root:

> f(A) = A^(-1/2)

from quantum chemistry. NTPoly is a massively parallel library that can be used
to compute a variety of matrix using polynomial expansions. Consider for example
the Taylor series expansion of a function *f(x)* .

> f(x) = f(0) + f'(0)x + f''(0)x^2/2! + ...

We can imagine expanding this from the function of a single variable, to a
function of a matrix:

> f(A) = f(0) + f'(0)A + f''(0)A^2/2! + ...

where matrices can be summed using matrix addition, and raised to a power
using matrix multiplication. At the heart of NTPoly are polynomial expansions
like this. We implement not only Taylor expansions, but also Chebyshev
polynomial expansions, and other specialized expansions based on the function
of interest.

When the input matrix *A* and the output matrix *f(A)* are sparse, we can
replace the dense matrix addition and multiplication routines with sparse
matrix routines. This allows us to use NTPoly to efficiently compute many
functions of sparse matrices.

Getting Start With Examples
--------------------------------------------------------------------------------
In the examples directory, there are a number of different example programs that
use NTPoly. You can check the ReadMe.md file in each example directory to
learn how to build and run each example. The simplest example is PremadeMatrix,
which includes sample output you can compare to.

Feature Outline
--------------------------------------------------------------------------------
The following features and methods have been implemented in NTPoly:

* General Polynomials
    * Standard Polynomials
    * Chebyshev Polynomials
    * Hermite Polynomials
* Transcendental Functions
    * Trigonometric Functions
    * Exponential and Logarithm
* Matrix Roots
    * Square Root and Inverse Square Root
    * Matrix *p* th Root
* Quantum Chemistry
    * Density Matrix Purification
    * Chemical Potential Calculation
    * Geometry Optimization
* Other
    * Matrix Inverse/Moore-Penrose Pseudo Inverse
    * Sign Function/Polar Decomposition
    * Load Balancing Matrices
    * File I/O

Citation
--------------------------------------------------------------------------------
A description of the techniques used in NTPoly can be found in the following
Computer Physics Communications paper:

> Dawson, William, and Takahito Nakajima. "Massively parallel sparse matrix
> function calculations with NTPoly." Computer Physics Communications (2017).

Please cite this paper in accordance to the practices in your field.

How To Contribute
--------------------------------------------------------------------------------
To begin contributing to NTPoly, take a look at the
[Wiki](https://github.com/william-dawson/NTPoly/wiki) pages. The
[Contributing Guide](https://github.com/william-dawson/NTPoly/blob/master/CONTRIBUTING.md)
provides an overview of best development practices. Additionally, there is a
[Adding New Functionality](https://github.com/william-dawson/NTPoly/wiki/Adding-New-Functionality-(Example))
page which documents how one would go about adding a matrix function to NTPoly.
