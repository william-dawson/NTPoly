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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  !! @param[out] MatV the eigenvalues.
  SUBROUTINE EigenDecomposition_ldr(MatA, MatV, MatW)
    !! Parameters
    TYPE(DMTYPE), INTENT(IN) :: MatA
    TYPE(DMTYPE), INTENT(INOUT) :: MatV
    TYPE(DMTYPE), INTENT(INOUT), OPTIONAL :: MatW
    !! Local variables
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
    DOUBLE PRECISION :: WORKTEMP
    INTEGER :: LWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: IWORKTEMP
    INTEGER :: LIWORK
    INTEGER :: INFO
    INTEGER :: II

    CALL ConstructEmptyMatrix(MatV,MatA%columns,MatA%rows)
    MatV%data = MatA%data

    N = SIZE(MatA%data,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL DSYEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, IWORKTEMP, &
         & LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP)
    ALLOCATE(WORK(LWORK))
    LIWORK = INT(IWORKTEMP)
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL DSYEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, IWORK, LIWORK, &
         & INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       CALL ConstructEmptyMatrix(MatW,1,MatA%rows)
       DO II = 1, N
          MatW%data(II,1) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)

  END SUBROUTINE EigenDecomposition_ldr

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  !! @param[out] MatV the eigenvalues.
  SUBROUTINE EigenDecomposition_ldc(MatA, MatV, MatW)
    !! Parameters
    TYPE(DMTYPE), INTENT(IN) :: MatA
    TYPE(DMTYPE), INTENT(INOUT) :: MatV
    TYPE(DMTYPE), INTENT(INOUT), OPTIONAL :: MatW
    !! Standard parameters
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: LWORK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
    INTEGER :: LRWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: LIWORK
    INTEGER :: INFO
    !! Temp
    COMPLEX*16 :: WORKTEMP
    DOUBLE PRECISION :: RWORKTEMP
    INTEGER :: IWORKTEMP
    INTEGER :: II

    CALL ConstructEmptyMatrix(MatV,MatA%columns,MatA%rows)
    MatV%data = MatA%data

    N = SIZE(MatA%data,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL ZHEEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, RWORKTEMP, &
         & LRWORK, IWORKTEMP, LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP)
    ALLOCATE(WORK(LWORK))
    LRWORK = INT(RWORKTEMP)
    ALLOCATE(RWORK(LRWORK))
    LIWORK = INT(IWORKTEMP)
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL ZHEEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, RWORK, LWORK, &
         & IWORK, LIWORK, INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       CALL ConstructEmptyMatrix(MatW,1,MatA%rows)
       DO II = 1, N
          MatW%data(II,1) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)
    DEALLOCATE(RWork)

  END SUBROUTINE EigenDecomposition_ldc

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
