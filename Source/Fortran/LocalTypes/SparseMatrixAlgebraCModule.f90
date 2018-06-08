!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE SparseMatrixAlgebraCModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE DenseMatrixCModule, ONLY : DenseMatrix_c, ConstructDenseFromSparse, &
       & ConstructSparseFromDense, MultiplyDense, DestructDenseMatrix
  USE MatrixMemoryPoolCModule, ONLY : MatrixMemoryPool_c, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE SparseMatrixCModule, ONLY: SparseMatrix_c, ConstructEmptySparseMatrix, &
       & DestructSparseMatrix, ConstructFromTripletList, CopySparseMatrix, &
       & TransposeSparseMatrix, PrintSparseMatrix
  USE SparseVectorCModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TripletListCModule, ONLY: TripletList_c, SortTripletList, &
       & ConstructTripletList, DestructTripletList

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE DenseMatrix_c
#define MPOOLTYPE MatrixMemoryPool_c
#define SMTYPE SparseMatrix_c
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c
#include "includes/SparseMatrixAlgebraImplementation.f90"

END MODULE SparseMatrixAlgebraCModule
