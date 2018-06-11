!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE MatrixDModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc, &
       & ConstructMatrixFromTripletList
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, AppendToTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldr
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldc
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE Matrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: ConstructMatrixDFromS
  PUBLIC :: ConstructMatrixSFromD
  PUBLIC :: CopyMatrix
  PUBLIC :: DestructMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SplitMatrix
  PUBLIC :: ComposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: EigenDecomposition
  PUBLIC :: MatrixNorm
  PUBLIC :: IncrementMatrix
  PUBLIC :: MultiplyMatrix
  PUBLIC :: TransposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrix_ldr
     MODULE PROCEDURE ConstructEmptyMatrix_ldc
  END INTERFACE
  INTERFACE ConstructMatrixDFromS
     MODULE PROCEDURE ConstructMatrixDFromS_ldr
     MODULE PROCEDURE ConstructMatrixDFromS_ldc
  END INTERFACE
  INTERFACE ConstructMatrixSFromD
     MODULE PROCEDURE ConstructMatrixSFromD_ldr
     MODULE PROCEDURE ConstructMatrixSFromD_ldc
  END INTERFACE
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_ldr
     MODULE PROCEDURE CopyMatrix_ldc
  END INTERFACE
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_ldr
     MODULE PROCEDURE DestructMatrix_ldc
  END INTERFACE
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_ldr
     MODULE PROCEDURE SplitMatrix_ldc
  END INTERFACE
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_ldr
     MODULE PROCEDURE ComposeMatrix_ldc
  END INTERFACE
  INTERFACE EigenDecomposition
     MODULE PROCEDURE EigenDecomposition_ldr
     MODULE PROCEDURE EigenDecomposition_ldc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ldr
     MODULE PROCEDURE MatrixNorm_ldc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ldr
     MODULE PROCEDURE IncrementMatrix_ldc
  END INTERFACE
  INTERFACE MultiplyMatrix
     MODULE PROCEDURE MultiplyMatrix_ldr
     MODULE PROCEDURE MultiplyMatrix_ldc
  END INTERFACE
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_ldr
     MODULE PROCEDURE TransposeMatrix_ldc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define DATATYPE REAL(NTREAL)
#define DMTYPE Matrix_ldr
#define SMTYPE Matrix_lsr
#define TTYPE Triplet_r
#define TLISTTYPE TripletList_r
#define ConstructEmptyMatrix ConstructEmptyMatrix_ldr
#define ConstructMatrixDFromS ConstructMatrixDFromS_ldr
#define ConstructMatrixSFromD ConstructMatrixSFromD_ldr
#define CopyMatrix CopyMatrix_ldr
#define DestructMatrix DestructMatrix_ldr
#define SplitMatrix SplitMatrix_ldr
#define ComposeMatrix ComposeMatrix_ldr
#define EigenDecomposition EigenDecomposition_ldr
#define MatrixNorm MatrixNorm_ldr
#define IncrementMatrix IncrementMatrix_ldr
#define MultiplyMatrix MultiplyMatrix_ldr
#define TransposeMatrix TransposeMatrix_ldr

#include "includes/MatrixDImpl.f90"

#undef ConstructEmptyMatrix
#undef ConstructMatrixDFromS
#undef ConstructMatrixSFromD
#undef CopyMatrix
#undef DestructMatrix
#undef SplitMatrix
#undef ComposeMatrix
#undef EigenDecomposition
#undef MatrixNorm
#undef IncrementMatrix
#undef MultiplyMatrix
#undef TransposeMatrix
#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

#define DATATYPE COMPLEX(NTCOMPLEX)
#define DMTYPE Matrix_ldc
#define SMTYPE Matrix_lsc
#define TTYPE Triplet_c
#define TLISTTYPE TripletList_c
#define ConstructEmptyMatrix ConstructEmptyMatrix_ldc
#define ConstructMatrixDFromS ConstructMatrixDFromS_ldc
#define ConstructMatrixSFromD ConstructMatrixSFromD_ldc
#define CopyMatrix CopyMatrix_ldc
#define DestructMatrix DestructMatrix_ldc
#define SplitMatrix SplitMatrix_ldc
#define ComposeMatrix ComposeMatrix_ldc
#define EigenDecomposition EigenDecomposition_ldc
#define MatrixNorm MatrixNorm_ldc
#define IncrementMatrix IncrementMatrix_ldc
#define MultiplyMatrix MultiplyMatrix_ldc
#define TransposeMatrix TransposeMatrix_ldc

#include "includes/MatrixDImpl.f90"

#undef ConstructEmptyMatrix
#undef ConstructMatrixDFromS
#undef ConstructMatrixSFromD
#undef CopyMatrix
#undef DestructMatrix
#undef SplitMatrix
#undef ComposeMatrix
#undef EigenDecomposition
#undef MatrixNorm
#undef IncrementMatrix
#undef MultiplyMatrix
#undef TransposeMatrix
#undef TLISTTYPE
#undef TTYPE
#undef SMTYPE
#undef DMTYPE
#undef DATATYPE

END MODULE MatrixDModule
