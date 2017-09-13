---
layout: page
title: About
permalink: /about/
---

## Basic Theory
The theory of matrix functions is a long studied branch of matrix algebra.
Matrix functions have a wide range of applications, including graph problems,
differential equations, and materials science. Common examples of matrix
functions include the matrix exponential:

$$ f(A) = e^A $$

from the study of networks, or the inverse square root:

$$ f(A) = A^{-\frac{1}{2}} $$

from quantum chemistry. NTPoly is a massively parallel library that can be used
to compute a variety of matrix using polynomial expansions. Consider for example
the Taylor series expansion of a function *f(x)* .

$$ f(x) = f(0) + f'(0)x + f''(0)\frac{x^2}{2!} + ...$$

We can imagine expanding this from the function of a single variable, to a
function of a matrix:

$$ f(A) = f(0) + f'(0)A + f''(0)\frac{A^2}{2!} + $$

where matrices can be summed using matrix addition, and raised to a power
using matrix multiplication. At the heart of NTPoly are polynomial expansions
like this. We implement not only Taylor expansions, but also Chebyshev
polynomial expansions, and other specialized expansions based on the function
of interest.

When the input matrix *A* and the output matrix *f(A)* are sparse, we can
replace the dense matrix addition and multiplication routines with sparse
matrix routines. This allows us to use NTPoly to efficiently compute many
functions of sparse matrices.

## How To Contribute

To begin contributing to NTPoly, take a look at the
[Wiki](https://github.com/william-dawson/NTPoly/wiki) pages. The
[Developer Guide](https://github.com/william-dawson/NTPoly/wiki/Developer-Guide)
provides an overview of best development practices. Additionally, there is a
[Adding New Functionality](https://github.com/william-dawson/NTPoly/wiki/Adding-New-Functionality-(Example))
page which documents how one would go about adding a matrix function to NTPoly.
