!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE MatrixSRAlgebraModule
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixDRModule, ONLY : Matrix_dr, ConstructDenseFromSparse, &
       & ConstructSparseFromDense, MultiplyDense, DestructDenseMatrix
  USE MatrixMemoryPoolRModule, ONLY : MatrixMemoryPool_r, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE MatrixSRModule, ONLY: Matrix_sr, ConstructEmptySparseMatrix, &
       & DestructSparseMatrix, ConstructFromTripletList, CopySparseMatrix, &
       & TransposeSparseMatrix, PrintSparseMatrix
  USE VectorSRModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TripletListRModule, ONLY: TripletList_r, SortTripletList, &
       & ConstructTripletList, DestructTripletList

#define DATATYPE REAL(NTREAL)
#define DMTYPE Matrix_dr
#define MPOOLTYPE MatrixMemoryPool_r
#define SMTYPE Matrix_sr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r

#include "includes/MatrixSAlgebraImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixSRAlgebraModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE MatrixSCAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixDCModule, ONLY : Matrix_dc, ConstructDenseFromSparse, &
       & ConstructSparseFromDense, MultiplyDense, DestructDenseMatrix
  USE MatrixMemoryPoolCModule, ONLY : MatrixMemoryPool_c, &
       & ConstructMatrixMemoryPool, DestructMatrixMemoryPool, &
       & CheckMemoryPoolValidity, SetPoolSparsity
  USE MatrixSCModule, ONLY: Matrix_sc, ConstructEmptySparseMatrix, &
       & DestructSparseMatrix, ConstructFromTripletList, CopySparseMatrix, &
       & TransposeSparseMatrix, PrintSparseMatrix
  USE VectorSCModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TripletListCModule, ONLY: TripletList_c, SortTripletList, &
       & ConstructTripletList, DestructTripletList

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE Matrix_dc
#define MPOOLTYPE MatrixMemoryPool_c
#define SMTYPE Matrix_sc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c

#include "includes/MatrixSAlgebraImpl.f90"

#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixSCAlgebraModule
