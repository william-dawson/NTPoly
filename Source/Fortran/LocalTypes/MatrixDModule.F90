!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE MatrixDModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE MatrixModule, ONLY : Matrix_l
  USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc
  USE TripletListModule, ONLY : TripletList_r, TripletList_c
  USE TripletModule, ONLY : Triplet_c, Triplet_r
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConvertDMatrixToS
  PUBLIC :: ConvertSMatrixToD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(Matrix_l), PUBLIC :: Matrix_ldr
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix.
  CONTAINS
    !! Construct/Destruct
    PROCEDURE :: InitEmpty => ConstructEmpty_ldr
    PROCEDURE :: Copy => Copy_ldr
    PROCEDURE :: Destruct => Destruct_ldr
    !! Basic Accessors
    PROCEDURE :: ExtractRow => ExtractRow_ldr
    PROCEDURE :: ExtractColumn => ExtractColumn_ldr
    !! ETC
    PROCEDURE :: Transpose => Transpose_ldr
    PROCEDURE :: Print => Print_ldr
    PROCEDURE :: PrintHeader => PrintHeader_ldr
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(Matrix_l), PUBLIC :: Matrix_ldc
    !> values of the matrix.
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: DATA
  CONTAINS
    !! Construct/Destruct
    PROCEDURE :: InitEmpty => ConstructEmpty_ldc
    PROCEDURE :: Copy => Copy_ldc
    PROCEDURE :: Destruct => Destruct_ldc
    !! Basic Accessors
    PROCEDURE :: ExtractRow => ExtractRow_ldc
    PROCEDURE :: ExtractColumn => ExtractColumn_ldc
    !! Algebra
    PROCEDURE :: Transpose => Transpose_ldc
    PROCEDURE :: Conjg => Conjg_ldc
    PROCEDURE :: Print => Print_ldc
    PROCEDURE :: PrintHeader => PrintHeader_ldc
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConvertDMatrixToS
     MODULE PROCEDURE ConvertDMatrixToS_r
     MODULE PROCEDURE ConvertDMatrixToS_c
  END INTERFACE
  INTERFACE ConvertSMatrixToD
     MODULE PROCEDURE ConvertSMatrixToD_r
     MODULE PROCEDURE ConvertSMatrixToD_c
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  !! @param[in] dense_matrix to convert.
  !! @param[out] sparse_matrix output matrix.
  !! @param[in] threshold_in value for pruning values to zero.
  PURE SUBROUTINE ConvertDMatrixToS_r(dense_matrix, sparse_matrix, &
       & threshold_in)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(IN) :: dense_matrix
    CLASS(Matrix_lsr), INTENT(INOUT) :: sparse_matrix
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_r) :: temporary
    TYPE(TripletList_r) :: temporary_list

    INCLUDE "dense_includes/ConstructMatrixSFromD.f90"

  END SUBROUTINE ConvertDMatrixToS_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  !! @param[in] dense_matrix to convert.
  !! @param[out] sparse_matrix output matrix.
  !! @param[in] threshold_in value for pruning values to zero.
  PURE SUBROUTINE ConvertDMatrixToS_c(dense_matrix, sparse_matrix, &
       & threshold_in)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(IN) :: dense_matrix
    CLASS(Matrix_lsc), INTENT(INOUT) :: sparse_matrix
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_c) :: temporary
    TYPE(TripletList_c) :: temporary_list

    INCLUDE "dense_includes/ConstructMatrixSFromD.f90"

  END SUBROUTINE ConvertDMatrixToS_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  !! @param[in] sparse_matrix a sparse matrix to convert.
  !! @param[inout] dense_matrix output. Must be preallocated.
  PURE SUBROUTINE ConvertSMatrixToD_r(sparse_matrix, dense_matrix)
    !! Parameters
    CLASS(Matrix_lsr), INTENT(IN) :: sparse_matrix
    CLASS(Matrix_ldr), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_r) :: temporary

    INCLUDE "dense_includes/ConstructMatrixDFromS.f90"

  END SUBROUTINE ConvertSMatrixToD_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  !! @param[in] sparse_matrix a sparse matrix to convert.
  !! @param[inout] dense_matrix output. Must be preallocated.
  PURE SUBROUTINE ConvertSMatrixToD_c(sparse_matrix, dense_matrix)
    !! Parameters
    CLASS(Matrix_lsc), INTENT(IN) :: sparse_matrix
    CLASS(Matrix_ldc), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_c) :: temporary

    INCLUDE "dense_includes/ConstructMatrixDFromS.f90"

  END SUBROUTINE ConvertSMatrixToD_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  !! @param[inout] this the matrix to construct.
  !! @param[in] columns of the matrix.
  !! @parma[in] rows of the matrix
  PURE SUBROUTINE ConstructEmpty_ldr(this, rows, columns, zero_in)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: rows
    INTEGER, INTENT(IN) :: columns
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "dense_includes/ConstructEmptyMatrix.f90"

  END SUBROUTINE ConstructEmpty_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  !! @param[inout] this the matrix to construct.
  !! @param[in] columns of the matrix.
  !! @parma[in] rows of the matrix
  PURE SUBROUTINE ConstructEmpty_ldc(this, rows, columns, zero_in)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: rows
    INTEGER, INTENT(IN) :: columns
    LOGICAL, INTENT(IN), OPTIONAL :: zero_in

    INCLUDE "dense_includes/ConstructEmptyMatrix.f90"

  END SUBROUTINE ConstructEmpty_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  !! @param[in] matA the matrix to copy.
  !! @param[inout] this = matA
  PURE SUBROUTINE Copy_ldr(this, matA)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(INOUT) :: this
    CLASS(Matrix_l), INTENT(IN) :: matA

    SELECT TYPE(matA)
    CLASS IS (Matrix_ldr)
       INCLUDE "dense_includes/CopyMatrix.f90"
    END SELECT

  END SUBROUTINE Copy_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  !! @param[in] matA the matrix to copy.
  !! @param[inout] this = matA
  PURE SUBROUTINE Copy_ldc(this, matA)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(INOUT) :: this
    CLASS(Matrix_l), INTENT(IN) :: matA

    SELECT TYPE(matA)
    CLASS IS (Matrix_ldc)
       INCLUDE "dense_includes/CopyMatrix.f90"
    END SELECT

  END SUBROUTINE Copy_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  !! @param[inout] this the matrix to delete.
  PURE SUBROUTINE Destruct_ldr(this)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(INOUT) :: this

    INCLUDE "dense_includes/DestructMatrix.f90"

  END SUBROUTINE Destruct_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  !! @param[inout] this the matrix to delete.
  PURE SUBROUTINE Destruct_ldc(this)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(INOUT) :: this

    INCLUDE "dense_includes/DestructMatrix.f90"

  END SUBROUTINE Destruct_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractRow_ldr(this, row_number, row_out)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    CLASS(Matrix_l), INTENT(INOUT) :: row_out

    SELECT TYPE(row_out)
    CLASS IS(Matrix_ldr)
      INCLUDE "dense_includes/ExtractMatrixRow.f90"
    END SELECT
  END SUBROUTINE ExtractRow_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractRow_ldc(this, row_number, row_out)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    CLASS(Matrix_l), INTENT(INOUT) :: row_out

    SELECT TYPE(row_out)
    CLASS IS(Matrix_ldc)
      INCLUDE "dense_includes/ExtractMatrixRow.f90"
    END SELECT
  END SUBROUTINE ExtractRow_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractColumn_ldr(this, column_number, column_out)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    CLASS(Matrix_l), INTENT(INOUT) :: column_out

    SELECT TYPE(column_out)
    CLASS IS(Matrix_ldr)
      INCLUDE "dense_includes/ExtractMatrixColumn.f90"
    END SELECT
  END SUBROUTINE ExtractColumn_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractColumn_ldc(this, column_number, column_out)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    CLASS(Matrix_l), INTENT(INOUT) :: column_out

    SELECT TYPE(column_out)
    CLASS IS(Matrix_ldc)
      INCLUDE "dense_includes/ExtractMatrixColumn.f90"
    END SELECT
  END SUBROUTINE ExtractColumn_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! @param[in] this = matA.T
  !! @param[out] matA the matrix to transpose.
  PURE SUBROUTINE Transpose_ldr(this, matA)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(INOUT)  :: this
    CLASS(Matrix_l), INTENT(IN) :: matA

    SELECT TYPE(matA)
    CLASS IS(Matrix_ldr)
       INCLUDE "dense_includes/TransposeMatrix.f90"
    END SELECT

  END SUBROUTINE Transpose_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! @param[in] this = matA.T
  !! @param[out] matA the matrix to transpose.
  PURE SUBROUTINE Transpose_ldc(this, matA)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(INOUT)  :: this
    CLASS(Matrix_l), INTENT(IN) :: matA

    SELECT TYPE(matA)
    CLASS IS(Matrix_ldc)
       INCLUDE "dense_includes/TransposeMatrix.f90"
    END SELECT

  END SUBROUTINE Transpose_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Every value in the matrix is changed into its complex conjugate.
  !! @param[inout] this the matrix to compute the complex conjugate of.
  PURE SUBROUTINE Conjg_ldc(this)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(INOUT) :: this

    this%data = CONJG(this%data)
  END SUBROUTINE Conjg_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE Print_ldr(this, file_name_in)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in

    INCLUDE "dense_includes/PrintMatrix.f90"

  END SUBROUTINE Print_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE Print_ldc(this, file_name_in)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in

    INCLUDE "dense_includes/PrintMatrix.f90"

  END SUBROUTINE Print_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the matrix market header for this matrix
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_handle an open file handle to write to.
  SUBROUTINE PrintHeader_ldr(this, file_handle)
    !! Parameters
    CLASS(Matrix_ldr), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,'(A)') "%%MatrixMarket matrix array real general"

  END SUBROUTINE PrintHeader_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out the matrix market header for this matrix
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_handle an open file handle to write to.
  SUBROUTINE PrintHeader_ldc(this, file_handle)
    !! Parameters
    CLASS(Matrix_ldc), INTENT(IN) :: this
    INTEGER, INTENT(INOUT) :: file_handle

    WRITE(file_handle,'(A)') "%%MatrixMarket matrix array complex general"

  END SUBROUTINE PrintHeader_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixDModule
! MODULE MatrixDModule
!   USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
!   USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc
!   USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
!        & AppendToTripletList, ConstructTripletList
!   USE TripletModule, ONLY : Triplet_r, Triplet_c
!   IMPLICIT NONE
!   PRIVATE
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A datatype for storing a dense matrix.
!   TYPE, PUBLIC :: Matrix_ldr
!      REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix.
!      INTEGER :: rows !< Matrix dimension: rows.
!      INTEGER :: columns !< Matrix dimension: columns.
!   END TYPE Matrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A datatype for storing a dense matrix.
!   TYPE, PUBLIC :: Matrix_ldc
!      !> values of the matrix.
!      COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: DATA
!      INTEGER :: rows !< Matrix dimension: rows.
!      INTEGER :: columns !< Matrix dimension: columns.
!    CONTAINS
!      FINAL :: DestructMatrix_ldc
!   END TYPE Matrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PUBLIC :: ConstructEmptyMatrix
!   PUBLIC :: ConstructMatrixDFromS
!   PUBLIC :: ConstructMatrixSFromD
!   PUBLIC :: CopyMatrix
!   PUBLIC :: DestructMatrix
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PUBLIC :: SplitMatrix
!   PUBLIC :: ComposeMatrix
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PUBLIC :: EigenDecomposition
!   PUBLIC :: MatrixNorm
!   PUBLIC :: IncrementMatrix
!   PUBLIC :: MultiplyMatrix
!   PUBLIC :: TransposeMatrix
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   INTERFACE Matrix_ldr
!      MODULE PROCEDURE ConstructEmptyMatrix_ldr
!   END INTERFACE
!   INTERFACE Matrix_ldc
!      MODULE PROCEDURE ConstructEmptyMatrix_ldc
!   END INTERFACE
!   INTERFACE ConstructEmptyMatrix
!      MODULE PROCEDURE ConstructEmptyMatrixSup_ldr
!      MODULE PROCEDURE ConstructEmptyMatrixSup_ldc
!   END INTERFACE
!   INTERFACE ConstructMatrixDFromS
!      MODULE PROCEDURE ConstructMatrixDFromS_ldr
!      MODULE PROCEDURE ConstructMatrixDFromS_ldc
!   END INTERFACE
!   INTERFACE ConstructMatrixSFromD
!      MODULE PROCEDURE ConstructMatrixSFromD_ldr
!      MODULE PROCEDURE ConstructMatrixSFromD_ldc
!   END INTERFACE
!   INTERFACE CopyMatrix
!      MODULE PROCEDURE CopyMatrix_ldr
!      MODULE PROCEDURE CopyMatrix_ldc
!   END INTERFACE
!   INTERFACE DestructMatrix
!      MODULE PROCEDURE DestructMatrix_ldr
!      MODULE PROCEDURE DestructMatrix_ldc
!   END INTERFACE
!   INTERFACE SplitMatrix
!      MODULE PROCEDURE SplitMatrix_ldr
!      MODULE PROCEDURE SplitMatrix_ldc
!   END INTERFACE
!   INTERFACE ComposeMatrix
!      MODULE PROCEDURE ComposeMatrix_ldr
!      MODULE PROCEDURE ComposeMatrix_ldc
!   END INTERFACE
!   INTERFACE EigenDecomposition
!      MODULE PROCEDURE EigenDecomposition_ldr
!      MODULE PROCEDURE EigenDecomposition_ldc
!   END INTERFACE
!   INTERFACE MatrixNorm
!      MODULE PROCEDURE MatrixNorm_ldr
!      MODULE PROCEDURE MatrixNorm_ldc
!   END INTERFACE
!   INTERFACE IncrementMatrix
!      MODULE PROCEDURE IncrementMatrix_ldr
!      MODULE PROCEDURE IncrementMatrix_ldc
!   END INTERFACE
!   INTERFACE MultiplyMatrix
!      MODULE PROCEDURE MultiplyMatrix_ldr
!      MODULE PROCEDURE MultiplyMatrix_ldc
!   END INTERFACE
!   INTERFACE TransposeMatrix
!      MODULE PROCEDURE TransposeMatrix_ldr
!      MODULE PROCEDURE TransposeMatrix_ldc
!   END INTERFACE
! CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE ConstructEmptyMatrixSup_ldr(this, rows, columns)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(INOUT) :: this
!     INTEGER, INTENT(IN) :: rows
!     INTEGER, INTENT(IN) :: columns
!
!     this = ConstructEmptyMatrix_ldr(rows, columns)
!   END SUBROUTINE ConstructEmptyMatrixSup_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Construct an empty dense matrix with a set number of rows and columns
!   !! @param[inout] this the matrix to construct.
!   !! @param[in] columns of the matrix.
!   !! @parma[in] rows of the matrix
!   PURE FUNCTION ConstructEmptyMatrix_ldr(rows, columns) RESULT(this)
!     !! Parameters
!     TYPE(Matrix_ldr) :: this
!     INTEGER, INTENT(IN) :: rows
!     INTEGER, INTENT(IN) :: columns
!
!     INCLUDE "dense_includes/ConstructEmptyMatrix.f90"
!
!   END FUNCTION ConstructEmptyMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A function that converts a sparse matrix to a dense matrix.
!   !! @param[in] sparse_matrix a sparse matrix to convert.
!   !! @param[inout] dense_matrix output. Must be preallocated.
!   PURE SUBROUTINE ConstructMatrixDFromS_ldr(sparse_matrix, dense_matrix)
!     !! Parameters
!     TYPE(Matrix_lsr), INTENT(IN) :: sparse_matrix
!     TYPE(Matrix_ldr), INTENT(INOUT) :: dense_matrix
!     !! Helper Variables
!     TYPE(Triplet_r) :: temporary
!
! #include "dense_includes/ConstructMatrixDFromS.f90"
!
!   END SUBROUTINE ConstructMatrixDFromS_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A function that converts a dense matrix to a sparse matrix.
!   !! @param[in] dense_matrix to convert.
!   !! @param[out] sparse_matrix output matrix.
!   !! @param[in] threshold_in value for pruning values to zero.
!   PURE SUBROUTINE ConstructMatrixSFromD_ldr(dense_matrix, sparse_matrix, &
!        & threshold_in)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: dense_matrix
!     TYPE(Matrix_lsr), INTENT(INOUT) :: sparse_matrix
!     REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
!     !! Local Variables
!     TYPE(Triplet_r) :: temporary
!     TYPE(TripletList_r) :: temporary_list
!
! #define SMTYPE Matrix_lsr
! #include "dense_includes/ConstructMatrixSFromD.f90"
! #undef SMTYPE
!
!   END SUBROUTINE ConstructMatrixSFromD_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Copy the matrix A into the B.
!   !! @param[in] matA the matrix to copy.
!   !! @param[inout] matB = matA
!   PURE SUBROUTINE CopyMatrix_ldr(matA, matB)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: matA
!     TYPE(Matrix_ldr), INTENT(INOUT) :: matB
!
! #include "dense_includes/CopyMatrix.f90"
!
!   END SUBROUTINE CopyMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Deallocate the memory associated with this matrix.
!   !! @param[inout] this the matrix to delete.
!   PURE SUBROUTINE DestructMatrix_ldr(this)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(INOUT) :: this
!
!     INCLUDE "dense_includes/DestructMatrix.f90"
!
!   END SUBROUTINE DestructMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> AXPY for dense matrices. B = B + alpha*A
!   !! @param[in] MatA is added
!   !! @param[inout] MatB is incremented.
!   !! @param[in] alpha_in a scaling parameter (optional).
!   PURE SUBROUTINE IncrementMatrix_ldr(MatA,MatB,alpha_in)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: MatA
!     TYPE(Matrix_ldr), INTENT(INOUT) :: MatB
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     !! Temporary
!     REAL(NTREAL) :: alpha
!
!     INCLUDE "dense_includes/IncrementMatrix.f90"
!
!   END SUBROUTINE IncrementMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the norm of a dense matrix.
!   !! Computes the Frobenius norm.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @result the norm of the matrix.
!   FUNCTION MatrixNorm_ldr(this) RESULT(norm)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: this
!     REAL(NTREAL) :: norm
!     !! Local Variables
!     CHARACTER, PARAMETER :: NORMC = 'F'
!     DOUBLE PRECISION, DIMENSION(this%rows) :: WORK
!     INTEGER :: M, N, LDA
!     !! Externel
!     DOUBLE PRECISION, EXTERNAL :: dlange
!
!     M = this%rows
!     N = this%columns
!     LDA = this%rows
!     norm = DLANGE(NORMC, M, N, this%data, LDA, WORK)
!   END FUNCTION MatrixNorm_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Transpose a dense matrix.
!   !! @param[in] matA the matrix to transpose.
!   !! @param[inout] matAT = matA^T.
!   PURE SUBROUTINE TransposeMatrix_ldr(matA, matAT)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: matA
!     TYPE(Matrix_ldr), INTENT(INOUT) :: matAT
!
! #include "dense_includes/TransposeMatrix.f90"
!
!   END SUBROUTINE TransposeMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Create a big matrix from an array of matrices by putting them one next
!   !! to another.
!   !! @param[in] mat_array 2d array of matrices to compose.
!   !! @param[in] block_rows the number of rows of the array of blocks.
!   !! @param[in] block_columns the number of columns of the array of blocks.
!   !! @param[out] out_matrix the composed matrix.
!   PURE SUBROUTINE ComposeMatrix_ldr(mat_array, block_rows, block_columns, &
!        & out_matrix)
!     !! Parameters
!     TYPE(Matrix_ldr), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
!          & mat_array
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     TYPE(Matrix_ldr), INTENT(INOUT) :: out_matrix
!
! #include "dense_includes/ComposeMatrix.f90"
!
!   END SUBROUTINE ComposeMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Split a sparse matrix into an array of sparse matrices.
!   !! @param[in] this the matrix to split.
!   !! @param[in] block_rows number of rows to split the matrix into.
!   !! @param[in] block_columns number of columns to split the matrix into.
!   !! @param[out] split_array a COLUMNxROW array for the output to go into.
!   !! @param[in] block_size_row_in specifies the size of the  rows.
!   !! @param[in] block_size_column_in specifies the size of the columns.
!   PURE SUBROUTINE SplitMatrix_ldr(this, block_rows, block_columns, &
!        & split_array, block_size_row_in, block_size_column_in)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: this
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     TYPE(Matrix_ldr), DIMENSION(:,:), INTENT(INOUT) :: split_array
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
!
! #include "dense_includes/SplitMatrix.f90"
!
!   END SUBROUTINE SplitMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A wrapper for multiplying two dense matrices.
!   !! @param[in] MatA the first matrix.
!   !! @param[in] MatB the second matrix.
!   !! @param[inout] MatC = MatA*MatB.
!   SUBROUTINE MultiplyMatrix_ldr(MatA,MatB,MatC)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: MatA
!     TYPE(Matrix_ldr), INTENT(IN) :: MatB
!     TYPE(Matrix_ldr), INTENT(INOUT) :: MatC
!     !! Local variables
!     CHARACTER, PARAMETER :: TRANSA = 'N'
!     CHARACTER, PARAMETER :: TRANSB = 'N'
!     INTEGER :: M
!     INTEGER :: N
!     INTEGER :: K
!     DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0
!     INTEGER :: LDA
!     INTEGER :: LDB
!     DOUBLE PRECISION, PARAMETER :: BETA = 0.0
!     INTEGER :: LDC
!
!     MatC = Matrix_ldr(MatA%rows,MatB%columns)
!
!     !! Setup Lapack
!     M = MatA%rows
!     N = MatB%columns
!     K = MatA%columns
!     LDA = M
!     LDB = K
!     LDC = M
!
!     CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
!          & LDB, BETA, MatC%data, LDC)
!
!   END SUBROUTINE MultiplyMatrix_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the eigenvectors of a dense matrix.
!   !! Wraps a standard dense linear algebra routine.
!   !! @param[in] MatA the matrix to decompose.
!   !! @param[out] MatV the eigenvectors.
!   !! @param[out] MatV the eigenvalues.
!   SUBROUTINE EigenDecomposition_ldr(MatA, MatV, MatW)
!     !! Parameters
!     TYPE(Matrix_ldr), INTENT(IN) :: MatA
!     TYPE(Matrix_ldr), INTENT(INOUT) :: MatV
!     TYPE(Matrix_ldr), INTENT(INOUT), OPTIONAL :: MatW
!     !! Local variables
!     CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
!     INTEGER :: N, LDA
!     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
!     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
!     DOUBLE PRECISION :: WORKTEMP
!     INTEGER :: LWORK
!     INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
!     INTEGER :: IWORKTEMP
!     INTEGER :: LIWORK
!     INTEGER :: INFO
!     INTEGER :: II
!
!     MatV = Matrix_ldr(MatA%rows,MatA%columns)
!     MatV%data = MatA%data
!
!     N = SIZE(MatA%data,DIM=1)
!     LDA = N
!
!     !! Allocations
!     ALLOCATE(W(N))
!
!     !! Determine the scratch space size
!     LWORK = -1
!     CALL DSYEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, IWORKTEMP, &
!          & LIWORK, INFO)
!     N = LDA
!     LWORK = INT(WORKTEMP)
!     ALLOCATE(WORK(LWORK))
!     LIWORK = INT(IWORKTEMP)
!     ALLOCATE(IWORK(LIWORK))
!
!     !! Run Lapack For Real
!     CALL DSYEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, IWORK, LIWORK, &
!          & INFO)
!
!     !! Extract Eigenvalues
!     IF (PRESENT(MatW)) THEN
!        MatW = Matrix_ldr(MatA%rows,1)
!        DO II = 1, N
!           MatW%data(II,1) = W(II)
!        END DO
!     END IF
!
!     !! Cleanup
!     DEALLOCATE(W)
!     DEALLOCATE(Work)
!
!   END SUBROUTINE EigenDecomposition_ldr
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PURE SUBROUTINE ConstructEmptyMatrixSup_ldc(this, rows, columns)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(INOUT) :: this
!     INTEGER, INTENT(IN) :: rows
!     INTEGER, INTENT(IN) :: columns
!
!     this = ConstructEmptyMatrix_ldc(rows, columns)
!   END SUBROUTINE ConstructEmptyMatrixSup_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Construct an empty dense matrix with a set number of rows and columns
!   !! @param[inout] this the matrix to construct.
!   !! @param[in] columns of the matrix.
!   !! @parma[in] rows of the matrix
!   PURE FUNCTION ConstructEmptyMatrix_ldc(rows, columns) RESULT(this)
!     !! Parameters
!     TYPE(Matrix_ldc) :: this
!     INTEGER, INTENT(IN) :: rows
!     INTEGER, INTENT(IN) :: columns
!
!     INCLUDE "dense_includes/ConstructEmptyMatrix.f90"
!
!   END FUNCTION ConstructEmptyMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A function that converts a sparse matrix to a dense matrix.
!   !! @param[in] sparse_matrix a sparse matrix to convert.
!   !! @param[inout] dense_matrix output. Must be preallocated.
!   PURE SUBROUTINE ConstructMatrixDFromS_ldc(sparse_matrix, dense_matrix)
!     !! Parameters
!     TYPE(Matrix_lsc), INTENT(IN) :: sparse_matrix
!     TYPE(Matrix_ldc), INTENT(INOUT) :: dense_matrix
!     !! Helper Variables
!     TYPE(Triplet_c) :: temporary
!
! #include "dense_includes/ConstructMatrixDFromS.f90"
!
!   END SUBROUTINE ConstructMatrixDFromS_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A function that converts a dense matrix to a sparse matrix.
!   !! @param[in] dense_matrix to convert.
!   !! @param[out] sparse_matrix output matrix.
!   !! @param[in] threshold_in value for pruning values to zero.
!   PURE SUBROUTINE ConstructMatrixSFromD_ldc(dense_matrix, sparse_matrix, &
!        & threshold_in)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: dense_matrix
!     TYPE(Matrix_lsc), INTENT(INOUT) :: sparse_matrix
!     REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
!     !! Local Variables
!     TYPE(Triplet_c) :: temporary
!     TYPE(TripletList_c) :: temporary_list
!
! #define SMTYPE Matrix_lsc
! #include "dense_includes/ConstructMatrixSFromD.f90"
! #undef SMTYPE
!
!   END SUBROUTINE ConstructMatrixSFromD_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Copy the matrix A into the B.
!   !! @param[in] matA the matrix to copy.
!   !! @param[inout] matB = matA
!   PURE SUBROUTINE CopyMatrix_ldc(matA, matB)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: matA
!     TYPE(Matrix_ldc), INTENT(INOUT) :: matB
!
! #include "dense_includes/CopyMatrix.f90"
!
!   END SUBROUTINE CopyMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Deallocate the memory associated with this matrix.
!   !! @param[inout] this the matrix to delete.
!   PURE SUBROUTINE DestructMatrix_ldc(this)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(INOUT) :: this
!
!     INCLUDE "dense_includes/DestructMatrix.f90"
!
!   END SUBROUTINE DestructMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> AXPY for dense matrices. B = B + alpha*A
!   !! @param[in] MatA is added
!   !! @param[inout] MatB is incremented.
!   !! @param[in] alpha_in a scaling parameter (optional).
!   PURE SUBROUTINE IncrementMatrix_ldc(MatA,MatB,alpha_in)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: MatA
!     TYPE(Matrix_ldc), INTENT(INOUT) :: MatB
!     REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
!     !! Temporary
!     REAL(NTREAL) :: alpha
!
!     INCLUDE "dense_includes/IncrementMatrix.f90"
!
!   END SUBROUTINE IncrementMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the norm of a dense matrix.
!   !! Computes the Frobenius norm.
!   !! @param[in] this the matrix to compute the norm of.
!   !! @result the norm of the matrix.
!   FUNCTION MatrixNorm_ldc(this) RESULT(norm)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: this
!     REAL(NTREAL) :: norm
!     !! Local Variables
!     CHARACTER, PARAMETER :: NORMC = 'F'
!     COMPLEX*16, DIMENSION(this%rows) :: WORK
!     INTEGER :: M, N, LDA
!     !! Externel
!     DOUBLE PRECISION, EXTERNAL :: zlange
!
!     M = this%rows
!     N = this%columns
!     LDA = this%rows
!     norm = ZLANGE(NORMC, M, N, this%data, LDA, WORK)
!   END FUNCTION MatrixNorm_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Transpose a dense matrix.
!   !! @param[in] matA the matrix to transpose.
!   !! @param[inout] matAT = matA^T.
!   PURE SUBROUTINE TransposeMatrix_ldc(matA, matAT)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: matA
!     TYPE(Matrix_ldc), INTENT(INOUT) :: matAT
!
! #include "dense_includes/TransposeMatrix.f90"
!
!   END SUBROUTINE TransposeMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Create a big matrix from an array of matrices by putting them one next
!   !! to another.
!   !! @param[in] mat_array 2d array of matrices to compose.
!   !! @param[in] block_rows the number of rows of the array of blocks.
!   !! @param[in] block_columns the number of columns of the array of blocks.
!   !! @param[out] out_matrix the composed matrix.
!   PURE SUBROUTINE ComposeMatrix_ldc(mat_array, block_rows, block_columns, &
!        & out_matrix)
!     !! Parameters
!     TYPE(Matrix_ldc), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
!          & mat_array
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     TYPE(Matrix_ldc), INTENT(INOUT) :: out_matrix
!
! #include "dense_includes/ComposeMatrix.f90"
!
!   END SUBROUTINE ComposeMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Split a sparse matrix into an array of sparse matrices.
!   !! @param[in] this the matrix to split.
!   !! @param[in] block_rows number of rows to split the matrix into.
!   !! @param[in] block_columns number of columns to split the matrix into.
!   !! @param[out] split_array a COLUMNxROW array for the output to go into.
!   !! @param[in] block_size_row_in specifies the size of the  rows.
!   !! @param[in] block_size_column_in specifies the size of the columns.
!   PURE SUBROUTINE SplitMatrix_ldc(this, block_rows, block_columns, &
!        & split_array, block_size_row_in, block_size_column_in)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: this
!     INTEGER, INTENT(IN) :: block_rows, block_columns
!     TYPE(Matrix_ldc), DIMENSION(:,:), INTENT(INOUT) :: split_array
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
!     INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
!
! #include "dense_includes/SplitMatrix.f90"
!
!   END SUBROUTINE SplitMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> A wrapper for multiplying two dense matrices.
!   !! @param[in] MatA the first matrix.
!   !! @param[in] MatB the second matrix.
!   !! @param[inout] MatC = MatA*MatB.
!   SUBROUTINE MultiplyMatrix_ldc(MatA,MatB,MatC)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: MatA
!     TYPE(Matrix_ldc), INTENT(IN) :: MatB
!     TYPE(Matrix_ldc), INTENT(INOUT) :: MatC
!     !! Local variables
!     CHARACTER, PARAMETER :: TRANSA = 'N'
!     CHARACTER, PARAMETER :: TRANSB = 'N'
!     INTEGER :: M
!     INTEGER :: N
!     INTEGER :: K
!     COMPLEX*16, PARAMETER :: ALPHA = 1.0
!     INTEGER :: LDA
!     INTEGER :: LDB
!     COMPLEX*16, PARAMETER :: BETA = 0.0
!     INTEGER :: LDC
!
!     MatC = Matrix_ldc(MatA%rows,MatB%columns)
!
!     !! Setup Lapack
!     M = MatA%rows
!     N = MatB%columns
!     K = MatA%columns
!     LDA = M
!     LDB = K
!     LDC = M
!
!     CALL ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
!          & LDB, BETA, MatC%data, LDC)
!
!   END SUBROUTINE MultiplyMatrix_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !> Compute the eigenvectors of a dense matrix.
!   !! Wraps a standard dense linear algebra routine.
!   !! @param[in] MatA the matrix to decompose.
!   !! @param[out] MatV the eigenvectors.
!   !! @param[out] MatV the eigenvalues.
!   SUBROUTINE EigenDecomposition_ldc(MatA, MatV, MatW)
!     !! Parameters
!     TYPE(Matrix_ldc), INTENT(IN) :: MatA
!     TYPE(Matrix_ldc), INTENT(INOUT) :: MatV
!     TYPE(Matrix_ldc), INTENT(INOUT), OPTIONAL :: MatW
!     !! Standard parameters
!     CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
!     INTEGER :: N, LDA
!     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
!     COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK
!     INTEGER :: LWORK
!     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
!     INTEGER :: LRWORK
!     INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
!     INTEGER :: LIWORK
!     INTEGER :: INFO
!     !! Temp
!     COMPLEX*16 :: WORKTEMP
!     DOUBLE PRECISION :: RWORKTEMP
!     INTEGER :: IWORKTEMP
!     INTEGER :: II
!
!     MatV = Matrix_ldc(MatA%rows,MatA%columns)
!     MatV%data = MatA%data
!
!     N = SIZE(MatA%data,DIM=1)
!     LDA = N
!
!     !! Allocations
!     ALLOCATE(W(N))
!
!     !! Determine the scratch space size
!     LWORK = -1
!     CALL ZHEEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, RWORKTEMP, &
!          & LRWORK, IWORKTEMP, LIWORK, INFO)
!     N = LDA
!     LWORK = INT(WORKTEMP)
!     ALLOCATE(WORK(LWORK))
!     LRWORK = INT(RWORKTEMP)
!     ALLOCATE(RWORK(LRWORK))
!     LIWORK = INT(IWORKTEMP)
!     ALLOCATE(IWORK(LIWORK))
!
!     !! Run Lapack For Real
!     CALL ZHEEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, RWORK, LRWORK, &
!          & IWORK, LIWORK, INFO)
!
!     !! Extract Eigenvalues
!     IF (PRESENT(MatW)) THEN
!        MatW = Matrix_ldc(MatA%rows, 1)
!        DO II = 1, N
!           MatW%data(II,1) = W(II)
!        END DO
!     END IF
!
!     !! Cleanup
!     DEALLOCATE(W)
!     DEALLOCATE(Work)
!     DEALLOCATE(RWork)
!
!   END SUBROUTINE EigenDecomposition_ldc
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END MODULE MatrixDModule
