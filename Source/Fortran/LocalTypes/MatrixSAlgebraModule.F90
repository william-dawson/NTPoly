!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for performing linear algebra using sparse matrices.
MODULE MatrixSAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixDModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructMatrixDFromS, &
       & ConstructMatrixSFromD, MultiplyMatrix, DestructMatrix
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
       & DestructMatrixMemoryPool, CheckMemoryPoolValidity, SetPoolSparsity
  USE MatrixSModule, ONLY: Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix, &
       & TransposeMatrix, PrintMatrix
  USE VectorSModule, ONLY : AddSparseVectors, DotSparseVectors, &
       & PairwiseMultiplyVectors
  USE TripletListModule, ONLY: TripletList_r, TripletList_c, SortTripletList, &
       & DestructTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ScaleMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: DotMatrix
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixColumnNorm
  PUBLIC :: MatrixNorm
  PUBLIC :: MatrixGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_lsr
     MODULE PROCEDURE ScaleMatrix_lsc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_lsr
     MODULE PROCEDURE IncrementMatrix_lsc
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_lsr
     MODULE PROCEDURE DotMatrix_lsc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsr
     MODULE PROCEDURE PairwiseMultiplyMatrix_lsc
  END INTERFACE
  INTERFACE MatrixMultiply
     MODULE PROCEDURE GemmMatrix_lsr
     MODULE PROCEDURE GemmMatrix_lsc
  END INTERFACE
  INTERFACE MatrixColumnNorm
     MODULE PROCEDURE MatrixColumnNorm_lsr
     MODULE PROCEDURE MatrixColumnNorm_lsc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_lsr
     MODULE PROCEDURE MatrixNorm_lsc
  END INTERFACE
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_lsr
     MODULE PROCEDURE MatrixGrandSum_lsc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define DATATYPE REAL(NTREAL)
#define DMTYPE Matrix_ldr
#define MPOOLTYPE MatrixMemoryPool_lr
#define SMTYPE Matrix_lsr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r
#define ScaleMatrix ScaleMatrix_lsr
#define IncrementMatrix IncrementMatrix_lsr
#define DotMatrix DotMatrix_lsr
#define PairwiseMultiplyMatrix PairwiseMultiplyMatrix_lsr
#define GemmMatrix GemmMatrix_lsr
#define MatrixColumnNorm MatrixColumnNorm_lsr
#define MatrixNorm MatrixNorm_lsr
#define MatrixGrandSum MatrixGrandSum_lsr
#define SparseBranch SparseBranch_lsr
#define DenseBranch DenseBranch_lsr
#define MultiplyBlock MultiplyBlock_lsr
#define PruneList PruneList_lsr

#include "includes/MatrixSAlgebraImpl.f90"

#undef ScaleMatrix
#undef IncrementMatrix
#undef DotMatrix
#undef PairwiseMultiplyMatrix
#undef GemmMatrix
#undef MatrixColumnNorm
#undef MatrixNorm
#undef MatrixGrandSum
#undef SparseBranch
#undef DenseBranch
#undef MultiplyBlock
#undef PruneList
#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE Matrix_ldc
#define MPOOLTYPE MatrixMemoryPool_lc
#define SMTYPE Matrix_lsc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c
#define ScaleMatrix ScaleMatrix_lsc
#define IncrementMatrix IncrementMatrix_lsc
#define DotMatrix DotMatrix_lsc
#define PairwiseMultiplyMatrix PairwiseMultiplyMatrix_lsc
#define GemmMatrix GemmMatrix_lsc
#define MatrixColumnNorm MatrixColumnNorm_lsc
#define MatrixNorm MatrixNorm_lsc
#define MatrixGrandSum MatrixGrandSum_lsc
#define SparseBranch SparseBranch_lsc
#define DenseBranch DenseBranch_lsc
#define MultiplyBlock MultiplyBlock_lsc
#define PruneList PruneList_lsc

#include "includes/MatrixSAlgebraImpl.f90"

#undef ScaleMatrix
#undef IncrementMatrix
#undef DotMatrix
#undef PairwiseMultiplyMatrix
#undef GemmMatrix
#undef MatrixColumnNorm
#undef MatrixNorm
#undef MatrixGrandSum
#undef SparseBranch
#undef DenseBranch
#undef MultiplyBlock
#undef PruneList
#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef MPOOLTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixSAlgebraModule
