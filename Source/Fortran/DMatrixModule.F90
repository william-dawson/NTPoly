!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DMatrixModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, &
       & ConstructMatrixFromTripletList
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & AppendToTripletList, ConstructTripletList, DestructTripletList
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldr
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix.
     INTEGER :: rows !< Matrix dimension: rows.
     INTEGER :: columns !< Matrix dimension: columns.
  END TYPE Matrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldc
     !> values of the matrix.
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: DATA
     INTEGER :: rows !< Matrix dimension: rows.
     INTEGER :: columns !< Matrix dimension: columns.
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
  PUBLIC :: MatrixNorm
  PUBLIC :: IncrementMatrix
  PUBLIC :: MultiplyMatrix
  PUBLIC :: TransposeMatrix
  PUBLIC :: EigenDecomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE Matrix_ldr
     MODULE PROCEDURE ConstructEmptyMatrix_ldr
  END INTERFACE Matrix_ldr
  INTERFACE Matrix_ldc
     MODULE PROCEDURE ConstructEmptyMatrix_ldc
  END INTERFACE Matrix_ldc
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrixSup_ldr
     MODULE PROCEDURE ConstructEmptyMatrixSup_ldc
  END INTERFACE ConstructEmptyMatrix
  INTERFACE ConstructMatrixDFromS
     MODULE PROCEDURE ConstructMatrixDFromS_ldr
     MODULE PROCEDURE ConstructMatrixDFromS_ldc
  END INTERFACE ConstructMatrixDFromS
  INTERFACE ConstructMatrixSFromD
     MODULE PROCEDURE ConstructMatrixSFromD_ldr
     MODULE PROCEDURE ConstructMatrixSFromD_ldc
  END INTERFACE ConstructMatrixSFromD
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_ldr
     MODULE PROCEDURE CopyMatrix_ldc
  END INTERFACE CopyMatrix
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_ldr
     MODULE PROCEDURE DestructMatrix_ldc
  END INTERFACE DestructMatrix
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_ldr
     MODULE PROCEDURE SplitMatrix_ldc
  END INTERFACE SplitMatrix
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_ldr
     MODULE PROCEDURE ComposeMatrix_ldc
  END INTERFACE ComposeMatrix
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ldr
     MODULE PROCEDURE MatrixNorm_ldc
  END INTERFACE MatrixNorm
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ldr
     MODULE PROCEDURE IncrementMatrix_ldc
  END INTERFACE IncrementMatrix
  INTERFACE MultiplyMatrix
     MODULE PROCEDURE MultiplyMatrix_ldr
     MODULE PROCEDURE MultiplyMatrix_ldc
  END INTERFACE MultiplyMatrix
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_ldr
     MODULE PROCEDURE TransposeMatrix_ldc
  END INTERFACE TransposeMatrix
  INTERFACE EigenDecomposition
     MODULE PROCEDURE EigenDecomposition_ldr
     MODULE PROCEDURE EigenDecomposition_ldc
  END INTERFACE EigenDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine wrapper for the empty constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSup_ldr(this, rows, columns)
    !> The matrix to construct
    TYPE(Matrix_ldr), INTENT(INOUT) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix
    INTEGER, INTENT(IN) :: columns

    CALL DestructMatrix(this)
    this = ConstructEmptyMatrix_ldr(rows, columns)
  END SUBROUTINE ConstructEmptyMatrixSup_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  PURE FUNCTION ConstructEmptyMatrix_ldr(rows, columns) RESULT(this)
    !> The matrix to construct.
    TYPE(Matrix_ldr) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix.
    INTEGER, INTENT(IN) :: columns

#include "dense_includes/ConstructEmptyMatrix.f90"

  END FUNCTION ConstructEmptyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  PURE SUBROUTINE ConstructMatrixDFromS_ldr(sparse_matrix, dense_matrix)
    !> The sparse matrix to convert.
    TYPE(Matrix_lsr), INTENT(IN) :: sparse_matrix
    !> Output. Must be preallocated.
    TYPE(Matrix_ldr), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_r) :: temporary

#include "dense_includes/ConstructMatrixDFromS.f90"

  END SUBROUTINE ConstructMatrixDFromS_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  PURE SUBROUTINE ConstructMatrixSFromD_ldr(dense_matrix, sparse_matrix, &
       & threshold_in)
    !> Matrix to convert.
    TYPE(Matrix_ldr), INTENT(IN) :: dense_matrix
    !> Output matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: sparse_matrix
    !> Value for pruning values to zero.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_r) :: temporary
    TYPE(TripletList_r) :: temporary_list

#include "dense_includes/ConstructMatrixSFromD.f90"

  END SUBROUTINE ConstructMatrixSFromD_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  PURE SUBROUTINE CopyMatrix_ldr(matA, matB)
    !> The matrix to copy.
    TYPE(Matrix_ldr), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_ldr), INTENT(INOUT) :: matB

#include "dense_includes/CopyMatrix.f90"

  END SUBROUTINE CopyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  PURE SUBROUTINE DestructMatrix_ldr(this)
    !> The matrix to delete.
    TYPE(Matrix_ldr), INTENT(INOUT) :: this

#include "dense_includes/DestructMatrix.f90"

  END SUBROUTINE DestructMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> AXPY for dense matrices. B = B + alpha*A
  PURE SUBROUTINE IncrementMatrix_ldr(MatA,MatB,alpha_in)
    !> MatA is added
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    !> MatB is incremented.
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatB
    !> A scaling parameter.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !! Temporary
    REAL(NTREAL) :: alpha

#include "dense_includes/IncrementMatrix.f90"

  END SUBROUTINE IncrementMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !> Computes the Frobenius norm.
  FUNCTION MatrixNorm_ldr(this) RESULT(norm)
    !! Parameters
    !> The matrix to compute the norm of.
    TYPE(Matrix_ldr), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    INTEGER :: II, JJ

    norm = 0
    DO II =1, this%rows
       DO JJ = 1,  this%columns
          norm = norm + this%DATA(II,JJ)**2
       END DO
    END DO
  END FUNCTION MatrixNorm_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a dense matrix.
  PURE SUBROUTINE TransposeMatrix_ldr(matA, matAT)
    !> matA the matrix to transpose.
    TYPE(Matrix_ldr), INTENT(IN) :: matA
    !> matAT = matA^T.
    TYPE(Matrix_ldr), INTENT(INOUT) :: matAT

#include "dense_includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_ldr(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_ldr), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
         & mat_array
    !> The composed matrix.
    TYPE(Matrix_ldr), INTENT(INOUT) :: out_matrix

#include "dense_includes/ComposeMatrix.f90"

  END SUBROUTINE ComposeMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_ldr(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_ldr), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A block_columns x block_rows array for the output to go into.
    TYPE(Matrix_ldr), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in

#include "dense_includes/SplitMatrix.f90"

  END SUBROUTINE SplitMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  SUBROUTINE MultiplyMatrix_ldr(MatA, MatB, MatC, IsATransposed_in, &
       & IsBTransposed_in)
    !> The first matrix.
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    !> The second matrix.
    TYPE(Matrix_ldr), INTENT(IN) :: MatB
    !> MatC = MatA*MatB.
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !! Local variables
    CHARACTER :: TRANSA
    CHARACTER :: TRANSB
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    DOUBLE PRECISION, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    !! Optional Parameters
    TRANSA = 'N'
    IF (PRESENT(IsATransposed_in)) THEN
       IF (IsATransposed_in) THEN
          TRANSA = 'T'
       END IF
    END IF
    TRANSB = 'N'
    IF (PRESENT(IsBTransposed_in)) THEN
       IF (IsBTransposed_in) THEN
          TRANSB = 'T'
       END IF
    END IF

    !! Setup Lapack
    IF (TRANSA .EQ. 'T') THEN
       M = MatA%columns
    ELSE
       M = MatA%rows
    END IF

    IF (TRANSB .EQ. 'T') THEN
       N = MatB%rows
    ELSE
       N = MatB%columns
    END IF

    IF (TRANSA .EQ. 'T') THEN
       K = MatA%rows
    ELSE
       K = MatA%columns
    END IF

    IF (TRANSA .EQ. 'T') THEN
       LDA = K
    ELSE
       LDA = M
    END IF

    IF (TRANSB .EQ. 'T') THEN
       LDB = N
    ELSE
       LDB = K
    END IF

    LDC = M

    !! Multiply
    CALL ConstructEmptyMatrix(MatC, M, N)
    CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%DATA, LDA, MatB%DATA, &
         & LDB, BETA, MatC%DATA, LDC)

  END SUBROUTINE MultiplyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !> Wraps a standard dense linear algebra routine.
  SUBROUTINE EigenDecomposition_ldr(MatA, MatV, MatW)
    !> MatA the matrix to decompose.
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    !> The eigenvectors.
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatV
    !> The eigenvalues.
    TYPE(Matrix_ldr), INTENT(INOUT), OPTIONAL :: MatW
    !! Local variables
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
    DOUBLE PRECISION, DIMENSION(1) :: WORKTEMP
    INTEGER :: LWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER, DIMENSION(1) :: IWORKTEMP
    INTEGER :: LIWORK
    INTEGER :: INFO
    INTEGER :: II

    MatV = Matrix_ldr(MatA%rows,MatA%columns)
    MatV%DATA = MatA%DATA

    N = SIZE(MatA%DATA,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL DSYEVD(JOB, UPLO, N, MatA%DATA, LDA, W, WORKTEMP, LWORK, IWORKTEMP, &
         & LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP(1))
    ALLOCATE(WORK(LWORK))
    LIWORK = INT(IWORKTEMP(1))
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL DSYEVD(JOB, UPLO, N, MatV%DATA, LDA, W, WORK, LWORK, IWORK, LIWORK, &
         & INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       MatW = Matrix_ldr(MatA%rows, MatA%columns)
       MatW%DATA = 0
       DO II = 1, N
          MatW%DATA(II,II) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)

  END SUBROUTINE EigenDecomposition_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine style wrapper for the constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSup_ldc(this, rows, columns)
    !> The matrix to construct.
    TYPE(Matrix_ldc), INTENT(INOUT) :: this
    !> The number of rows of the matrix.
    INTEGER, INTENT(IN) :: rows
    !> The number of columns o the matrix.
    INTEGER, INTENT(IN) :: columns

    CALL DestructMatrix(this)
    this = ConstructEmptyMatrix_ldc(rows, columns)
  END SUBROUTINE ConstructEmptyMatrixSup_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  PURE FUNCTION ConstructEmptyMatrix_ldc(rows, columns) RESULT(this)
    !> The matrix to construct.
    TYPE(Matrix_ldc) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix.
    INTEGER, INTENT(IN) :: columns

#include "dense_includes/ConstructEmptyMatrix.f90"

  END FUNCTION ConstructEmptyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  PURE SUBROUTINE ConstructMatrixDFromS_ldc(sparse_matrix, dense_matrix)
    !> The sparse matrix to convert.
    TYPE(Matrix_lsc), INTENT(IN) :: sparse_matrix
    !> Dense matrix output. Must be preallocated.
    TYPE(Matrix_ldc), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_c) :: temporary

#include "dense_includes/ConstructMatrixDFromS.f90"

  END SUBROUTINE ConstructMatrixDFromS_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  PURE SUBROUTINE ConstructMatrixSFromD_ldc(dense_matrix, sparse_matrix, &
       & threshold_in)
    !> The matrix to convert.
    TYPE(Matrix_ldc), INTENT(IN) :: dense_matrix
    !> The sparse output matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: sparse_matrix
    !> Value for pruning values to zero.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_c) :: temporary
    TYPE(TripletList_c) :: temporary_list

#include "dense_includes/ConstructMatrixSFromD.f90"

  END SUBROUTINE ConstructMatrixSFromD_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  PURE SUBROUTINE CopyMatrix_ldc(matA, matB)
    !> The matrix to copy.
    TYPE(Matrix_ldc), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_ldc), INTENT(INOUT) :: matB

#include "dense_includes/CopyMatrix.f90"

  END SUBROUTINE CopyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  PURE SUBROUTINE DestructMatrix_ldc(this)
    !> This the matrix to delete.
    TYPE(Matrix_ldc), INTENT(INOUT) :: this

#include "dense_includes/DestructMatrix.f90"

  END SUBROUTINE DestructMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> AXPY for dense matrices. B = B + alpha*A
  PURE SUBROUTINE IncrementMatrix_ldc(MatA,MatB,alpha_in)
    !> MatA is added
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    !> MatB is incremented.
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatB
    !> A scaling parameter.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !! Temporary
    REAL(NTREAL) :: alpha

#include "dense_includes/IncrementMatrix.f90"

  END SUBROUTINE IncrementMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !> Computes the Frobenius norm.
  FUNCTION MatrixNorm_ldc(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_ldc), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    INTEGER :: II, JJ
    COMPLEX(NTCOMPLEX) :: val, conjval

    norm = 0
    DO II =1, this%rows
       DO JJ = 1,  this%columns
          val = this%DATA(II,JJ)
          conjval = CONJG(val)
          norm = norm + REAL(val*conjval,KIND=NTREAL)
       END DO
    END DO
  END FUNCTION MatrixNorm_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a dense matrix.
  PURE SUBROUTINE TransposeMatrix_ldc(matA, matAT)
    !> The matrix to transpose.
    TYPE(Matrix_ldc), INTENT(IN) :: matA
    !> matAT = matA^T.
    TYPE(Matrix_ldc), INTENT(INOUT) :: matAT

#include "dense_includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_ldc(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_ldc), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
         & mat_array
    !> The composed matrix.
    TYPE(Matrix_ldc), INTENT(INOUT) :: out_matrix

#include "dense_includes/ComposeMatrix.f90"

  END SUBROUTINE ComposeMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_ldc(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_ldc), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A COLUMNxROW array for the output to go into.
    TYPE(Matrix_ldc), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in

#include "dense_includes/SplitMatrix.f90"

  END SUBROUTINE SplitMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  SUBROUTINE MultiplyMatrix_ldc(MatA, MatB, MatC, IsATransposed_in, &
       & IsBTransposed_in)
    !> The first matrix.
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    !> The second matrix.
    TYPE(Matrix_ldc), INTENT(IN) :: MatB
    !> MatC = MatA*MatB.
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatC
    !> True if A is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsATransposed_in
    !> True if B is already transposed.
    LOGICAL, OPTIONAL, INTENT(IN) :: IsBTransposed_in
    !! Local variables
    CHARACTER :: TRANSA
    CHARACTER :: TRANSB
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    COMPLEX*16, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    COMPLEX*16, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    !! Optional Parameters
    TRANSA = 'N'
    IF (PRESENT(IsATransposed_in)) THEN
       IF (IsATransposed_in) THEN
          TRANSA = 'T'
       END IF
    END IF
    TRANSB = 'N'
    IF (PRESENT(IsBTransposed_in)) THEN
       IF (IsBTransposed_in) THEN
          TRANSB = 'T'
       END IF
    END IF

    !! Setup Lapack
    IF (TRANSA .EQ. 'T') THEN
       M = MatA%columns
    ELSE
       M = MatA%rows
    END IF

    IF (TRANSB .EQ. 'T') THEN
       N = MatB%rows
    ELSE
       N = MatB%columns
    END IF

    IF (TRANSA .EQ. 'T') THEN
       K = MatA%rows
    ELSE
       K = MatA%columns
    END IF

    IF (TRANSA .EQ. 'T') THEN
       LDA = K
    ELSE
       LDA = M
    END IF

    IF (TRANSB .EQ. 'T') THEN
       LDB = N
    ELSE
       LDB = K
    END IF

    LDC = M

    !! Multiply
    CALL ConstructEmptyMatrix(MatC, M, N)
    CALL ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%DATA, LDA, MatB%DATA, &
         & LDB, BETA, MatC%DATA, LDC)

  END SUBROUTINE MultiplyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !> Wraps a standard dense linear algebra routine.
  SUBROUTINE EigenDecomposition_ldc(MatA, MatV, MatW)
    !> The matrix to decompose.
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    !> The eigenvectors.
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatV
    !> The eigenvalues.
    TYPE(Matrix_ldc), INTENT(INOUT), OPTIONAL :: MatW
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
    COMPLEX*16, DIMENSION(1) :: WORKTEMP
    DOUBLE PRECISION, DIMENSION(1) :: RWORKTEMP
    INTEGER, DIMENSION(1) :: IWORKTEMP
    INTEGER :: II

    MatV = Matrix_ldc(MatA%rows,MatA%columns)
    MatV%DATA = MatA%DATA

    N = SIZE(MatA%DATA,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL ZHEEVD(JOB, UPLO, N, MatA%DATA, LDA, W, WORKTEMP, LWORK, RWORKTEMP, &
         & LRWORK, IWORKTEMP, LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP(1))
    ALLOCATE(WORK(LWORK))
    LRWORK = INT(RWORKTEMP(1))
    ALLOCATE(RWORK(LRWORK))
    LIWORK = INT(IWORKTEMP(1))
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL ZHEEVD(JOB, UPLO, N, MatV%DATA, LDA, W, WORK, LWORK, RWORK, LRWORK, &
         & IWORK, LIWORK, INFO)

    !! Extract Eigenvalues
    IF (PRESENT(MatW)) THEN
       MatW = Matrix_ldc(MatA%rows, MatA%columns)
       MatW%DATA = 0
       DO II = 1, N
          MatW%DATA(II,II) = W(II)
       END DO
    END IF

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)
    DEALLOCATE(RWork)

  END SUBROUTINE EigenDecomposition_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DMatrixModule
