---
layout: post
title:  "Added Some Features That Were Basically Always There"
date:   2018-02-08 11:00:00 +0900
categories: blog update
---

I've added two new matrix functions to NTPoly. First, the 
[Moore-Penrose pseudoinverse](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse)
which can be used to approximate the inverse of a singular matrix. The second
is the
[Polar Decomposition](https://en.wikipedia.org/wiki/Polar_decomposition) of
a matrix.

Both of these functions can be implemented using already existing routines in
NTPoly. The Moore-Penrose pseudoinverse can be computed using Hotelling's
method, we just change the criteria for convergence. The following paper
by Li Weiguo gives a nice overview:

> Weiguo, Li, Li Juan, and Qiao Tiantian. "A family of iterative methods 
> for computing Moore–Penrose inverse of a matrix." Linear Algebra and 
> Its Applications 438.1 (2013): 47-56.

The Polar Decomposition can be computed using a Newton-Schultz iteration, 
identical to the matrix sign function. However, in this case it is necessary 
to transpose the matrix in the iteration. The Polar Decomposition is
closely related to the matrix sign function, so I've bundled them into the
same module.

> Higham, Nicholas J. "The matrix sign decomposition and its relation to 
> the polar decomposition." Linear Algebra and its Applications 212 (1994): 
> 3-20.

I made these simple modifications for the latest version of NTPoly 
available on Github.

