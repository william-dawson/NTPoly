!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE MatrixMarketModule, ONLY : MM_REAL, MM_COMPLEX
  USE TripletListModule, ONLY : TripletList
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for all kinds of local matrices.
  TYPE, ABSTRACT, PUBLIC :: Matrix_l
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
   CONTAINS
     !! Construct/Destruct
     PROCEDURE(ConstructEmptyMatrix_l), DEFERRED :: InitEmpty
     PROCEDURE(ConstructFromFile_l), DEFERRED :: InitFromFile
     PROCEDURE(DestructMatrix_l), DEFERRED :: Destruct
     PROCEDURE(CopyMatrix_l), DEFERRED :: Copy
     !! Basic Accessors
     PROCEDURE :: GetRows => GetRows_l
     PROCEDURE :: GetColumns => GetColumns_l
     PROCEDURE(ExtractMatrixRow_l), DEFERRED :: ExtractRow
     PROCEDURE(ExtractMatrixColumn_l), DEFERRED :: ExtractColumn
     !! Algebra
     PROCEDURE(ScaleMatrix_l), DEFERRED :: Scale
     PROCEDURE(IncrementMatrix_l), DEFERRED :: Increment
     PROCEDURE(PairwiseMultiplyMatrix_l), DEFERRED :: PairwiseMultiply
     PROCEDURE(MatrixColumnNorm_l), DEFERRED :: ColumnNorm
     PROCEDURE(MatrixNorm_l), DEFERRED :: Norm
     !! ETC
     PROCEDURE(ConvertToTripletList_l), DEFERRED :: ConvertToTripletList
     PROCEDURE(TransposeMatrix_l), DEFERRED :: Transpose
     PROCEDURE(PrintMatrix_l), DEFERRED :: Print
     PROCEDURE(PrintMatrixHeader_l), DEFERRED :: PrintHeader
  END TYPE Matrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
     PURE SUBROUTINE ConstructEmptyMatrix_l(this, rows, columns, zero_in)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       CLASS(Matrix_l), INTENT(INOUT) :: this
       INTEGER, INTENT(IN) :: columns, rows
       LOGICAL, INTENT(IN), OPTIONAL :: zero_in
     END SUBROUTINE ConstructEmptyMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE ConstructFromFile_l(this, file_name)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT) :: this
       CHARACTER(len=*), INTENT(IN) :: file_name
     END SUBROUTINE ConstructFromFile_l
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Explicitly destruct a sparse matrix.
     !! @param[inout] this the matrix to free up
     PURE SUBROUTINE DestructMatrix_l(this)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT) :: this
     END SUBROUTINE DestructMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Copy a sparse matrix in a safe way.
     !! @param[in] matA matrix to copy
     !! @param[inout] this = matA
     PURE SUBROUTINE CopyMatrix_l(this, matA)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT) :: this
       CLASS(Matrix_l), INTENT(IN) :: matA
     END SUBROUTINE CopyMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Extract a row from the matrix into the compressed vector representation.
     !! @param[in] this the matrix to extract from.
     !! @param[in] row_number the row to extract
     !! @param[out] row_out the matrix representing that row
     PURE SUBROUTINE ExtractMatrixRow_l(this, row_number, row_out)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: row_number
       CLASS(Matrix_l), INTENT(INOUT) :: row_out
     END SUBROUTINE ExtractMatrixRow_l
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Extract a column from the matrix into the compressed vector representation.
     !! @param[in] this the matrix to extract from.
     !! @param[in] column_number the row to extract
     !! @param[out] column_out the matrix representing that row
     PURE SUBROUTINE ExtractMatrixColumn_l(this, column_number, column_out)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: column_number
       CLASS(Matrix_l), INTENT(INOUT) :: column_out
     END SUBROUTINE ExtractMatrixColumn_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Transpose a sparse matrix and return it in a separate matrix.
     !! @param[in] this = matA.T
     !! @param[out] matA the matrix to transpose.
     PURE SUBROUTINE TransposeMatrix_l(this, matA)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT)  :: this
       CLASS(Matrix_l), INTENT(IN) :: matA
     END SUBROUTINE TransposeMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Construct a triplet list from a matrix.
     !! @param[in] this the matrix to construct the triplet list from.
     !! @param[out] triplet_list the triplet list we created.
     PURE SUBROUTINE ConvertToTripletList_l(this, triplet_list)
       USE TripletListModule, ONLY : TripletList
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       CLASS(TripletList), INTENT(INOUT) :: Triplet_list
     END SUBROUTINE ConvertToTripletList_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Print out a sparse matrix.
     !! @param[in] this the matrix to be printed.
     !! @param[in] file_name_in optionally, you can pass a file to print to.
     SUBROUTINE PrintMatrix_l(this, file_name_in)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
     END SUBROUTINE PrintMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Print out a sparse matrix.
     !! @param[in] this the matrix to be printed.
     !! @param[in] file_name_in optionally, you can pass a file to print to.
     SUBROUTINE PrintMatrixHeader_l(this, file_handle)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER, INTENT(INOUT) :: file_handle
     END SUBROUTINE PrintMatrixHeader_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Will scale a sparse matrix by a constant.
     !! @param[inout] matA Matrix A.
     !! @param[in] constant scale factor.
     PURE SUBROUTINE ScaleMatrix_l(matA,constant)
       USE DataTypesModule, ONLY : NTREAL
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT) :: matA
       REAL(NTREAL), INTENT(IN) :: constant
     END SUBROUTINE ScaleMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> this = alpha*Matrix A + this (AXPY).
     !! This will utilize the sparse vector addition routine.
     !! @param[in,out] this.
     !! @param[in] matA Matrix A.
     !! @param[in] alpha_in multiplier (optional, default=1.0)
     !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
     PURE SUBROUTINE IncrementMatrix_l(this, matA, alpha_in, threshold_in)
       USE DataTypesModule, ONLY : NTREAL
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT)  :: this
       CLASS(Matrix_l), INTENT(IN) :: matA
       REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
       REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
     END SUBROUTINE IncrementMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Pairwise Multiply two matrices.
     !! This will utilize the sparse vector pairwise routine.
     !! @param[in] matA Matrix A.
     !! @param[in] matB Matrix B.
     !! @param[in,out] this = MatA mult MatB.
     PURE SUBROUTINE PairwiseMultiplyMatrix_l(this, matA, matB)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT)  :: this
       CLASS(Matrix_l), INTENT(IN) :: matA
       CLASS(Matrix_l), INTENT(IN) :: matB
     END SUBROUTINE PairwiseMultiplyMatrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Compute the norm of a sparse matrix along the columns.
     !! @param[in] this the matrix to compute the norm of.
     !! @param[out] norm_per_column the norm value for each column in this matrix.
     PURE SUBROUTINE MatrixColumnNorm_l(this, norm_per_column)
       USE DataTypesModule, ONLY : NTREAL
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: norm_per_column
     END SUBROUTINE MatrixColumnNorm_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Compute the 1 norm of a sparse matrix.
     !! @param[in] this the matrix to compute the norm of.
     !! @result norm the matrix.
     PURE FUNCTION MatrixNorm_l(this) RESULT(norm)
       USE DataTypesModule, ONLY : NTREAL
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       REAL(NTREAL) :: norm
     END FUNCTION MatrixNorm_l
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetRows_l(this) RESULT(rows)
    !! Parameters
    CLASS(Matrix_l), INTENT(IN) :: this
    INTEGER :: rows

    rows = this%rows
  END FUNCTION GetRows_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetColumns_l(this) RESULT(columns)
    !! Parameters
    CLASS(Matrix_l), INTENT(IN) :: this
    INTEGER :: columns

    columns = this%columns
  END FUNCTION GetColumns_l
END MODULE MatrixModule
