!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE SparseMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL
  USE DenseMatrixModule, ONLY : DenseMatrix_t, ConstructDenseFromSparse, &
       & ConstructSparseFromDense, MultiplyDense, DestructDenseMatrix
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_t, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE SparseMatrixModule, ONLY: SparseMatrix_t, ConstructEmptySparseMatrix, &
       & DestructSparseMatrix, ConstructFromTripletList, CopySparseMatrix, &
       & TransposeSparseMatrix, PrintSparseMatrix
  USE SparseVectorModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TripletListModule, ONLY: TripletList_t, SortTripletList, &
       & ConstructTripletList, DestructTripletList

#define DATATYPE REAL(NTREAL)
#define DMTYPE DenseMatrix_t
#define MPOOLTYPE MatrixMemoryPool_t
#define SMTYPE SparseMatrix_t
#define TTYPE Triplet_t
#define TLISTTYPE TripletList_t
#include "includes/SparseMatrixAlgebraImplementation.f90"

END MODULE SparseMatrixAlgebraModule
