# Overlap Matrix Example
In this example, we will show an example of inverting an overlap matrix in
NTPoly. The main point of this example is to show you how to divide up
a matrix.

## Background

Imagine we are writing a program that computes an overlap matrix, and then
uses NTPoly to compute the inverse square root. Let's begin by imagining
we have 6 MPI processes, arranged in a 2x3 grid.

|        | Row 1      | Row 2     |
|--------| -----------| --------- |
Column 1 | Process 1  | Process 2 |
Column 2 | Process 3  | Process 4 |
Column 3 | Process 5  | Process 6 |

We can also divide up the overlap matrix in exactly the same way. So
the 6x6 matrix below:

|                         |
|-------------------------|
| a11 a12 a13 a14 a15 a16 |
| a21 a22 a23 a24 a25 a26 |
| a31 a32 a33 a34 a35 a36 |
| a41 a42 a43 a44 a45 a46 |
| a51 a52 a53 a54 a55 a56 |
| a61 a62 a63 a64 a65 a66 |

is divided up like this:

|          | Row 1                   | Row 2                   |
|----------|-------------------------|-------------------------|
| Column 1 | a11 a12 a13 a21 a22 a23 | a14 a15 a16 a24 a25 a26 |
| Column 2 | a31 a32 a33 a41 a42 a43 | a34 a35 a36 a44 a45 a46 |
| Column 3 | a51 a52 a53 a61 a62 a63 | a54 a55 a56 a64 a65 a66 |

Then process 1 is responsible for computing matrix elements:

| Process 1               |
|-------------------------|
| a11 a12 a13 a21 a22 a23 |

Once those elements are computed, we can construct the distributed matrix, and
then call the solver.

## Code Outline

The following steps are carried out.
1. Process the input parameters.
2. Construct the process grid.
3. Use the process grid to split up the matrix.
4. Compute the local elements of the overlap matrix.
5. Fill in the overlap matrix.
6. Solver for the inverse square root.
7. Write to file.

## Build System

Build with
mpif90 main.f90 -o example \
  -I../../Build/include \
  -L../../Build/lib -lNTPoly -fopenmp -llapack

And then run with:
mpirun -np 1 ./example \
--process_rows 1 --process_columns 1 --process_slices 1 \
--threshold 1e-6 --convergence_threshold 1e-5 \
--basis_functions 100

## Input Parameters

The input parameters are:
- process_rows how many process rows to use
- process_columns how many process columns to use
- process_slices for now leave this as 1.
- threshold any value less than this are set to zero.
- convergence_threshold determines when we consider the solver finished.
- basis_functions how big the overlap matrix is.

## Construct The Process Grid

Once all the input is read, it is mandatory to call ConstructProcessGrid. After
that it is safe to use NTPoly.

## Divide Up The Matrix

Next we need to divide up the matrix. First, we compute how many matrix rows
and columns are stored locally. This is just the total matrix size divided
by the number of process rows/columns.

Second, we compute the starting and ending rows/columns. In the ProcessGridModule,
there are global variables "my_row" and "my_column" which gives the xy coordinates
of a given process. This information can be used to compute the starting row
and column, and then we just add to those values the number of local rows/columns.

## Compute The Local Matrix Elements

Now that we know the local rows and columns of the matrix, we can compute
the local matrix elements. We do that by looping over the local rows and columns,
and calling the "ComputeIntegral" function. This ComputeIntegral function is of
course just a dummy function.

Since we don't know a priori which integrals can be neglected, what we do is
compute each one, and then test if the value is bigger than some threshold.
If it is, we add it to a list.  Finally, we convert that list to an NTPoly
TripletList type, and use it to fill the distributed matrix.

## The Rest

The program ends by following the same steps as in the other examples. The
solver parameters are set up, a call to a solver is made, and the output
matrix is printed to file. We also print the overlap matrix for reference.
These matrices can be visualized with the included python script.
