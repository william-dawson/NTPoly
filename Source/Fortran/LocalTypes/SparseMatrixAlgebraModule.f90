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

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE SparseMatrixAlgebraModule

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

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE SparseMatrixAlgebraCModule
