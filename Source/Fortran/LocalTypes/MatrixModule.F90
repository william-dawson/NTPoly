!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE MatrixModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE MatrixMarketModule, ONLY : MM_REAL, MM_COMPLEX
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
     PROCEDURE(ConstructEmptyMatrix_l), DEFERRED :: ConstructEmptyMatrix
     PROCEDURE(DestructMatrix_l), DEFERRED :: DestructMatrix
     PROCEDURE(CopyMatrix_l), DEFERRED :: CopyMatrix
     !! Basic Accessors
     PROCEDURE(GetMatrixRows_l), DEFERRED :: GetMatrixRows
     PROCEDURE(GetMatrixColumns_l), DEFERRED :: GetMatrixColumns
     PROCEDURE(ExtractMatrixRow_l), DEFERRED :: ExtractMatrixRow
     PROCEDURE(ExtractMatrixColumn_l), DEFERRED :: ExtractMatrixColumn
     !! Algebra
     PROCEDURE(ScaleMatrix_l), DEFERRED :: ScaleMatrix
     PROCEDURE(IncrementMatrix_l), DEFERRED :: IncrementMatrix
     PROCEDURE(PairwiseMultiplyMatrix_l), DEFERRED :: PairwiseMultiplyMatrix
     PROCEDURE(MatrixColumnNorm_l), DEFERRED :: MatrixColumnNorm
     PROCEDURE(MatrixNorm_l), DEFERRED :: MatrixNorm
     !! ETC
     PROCEDURE(TransposeMatrix_l), DEFERRED :: TransposeMatrix
     PROCEDURE(PrintMatrix_l), DEFERRED :: PrintMatrix
  END TYPE Matrix_l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
     PURE SUBROUTINE ConstructEmptyMatrix_l(this, rows, columns, zero_in)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: columns, rows
       LOGICAL, INTENT(IN), OPTIONAL :: zero_in
     END SUBROUTINE ConstructEmptyMatrix_l
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE ConstructMatrixFromFile_l(this, file_name)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(INOUT) :: this
       CHARACTER(len=*), INTENT(IN) :: file_name
     END SUBROUTINE ConstructMatrixFromFile_l
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
     !! @param[inout] matB = matA
     PURE SUBROUTINE CopyMatrix_l(matA, matB)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: matA
       CLASS(Matrix_l), INTENT(INOUT) :: matB
     END SUBROUTINE CopyMatrix_l
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Get the number of rows of a matrix.
     !! @param[in] this the matrix.
     !! @result number of rows.
     PURE FUNCTION GetMatrixRows_l(this) RESULT(rows)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER :: rows
     END FUNCTION GetMatrixRows_l
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Get the number of columns of a matrix.
     !! @param[in] this the matrix.
     !! @result number of columns.
     PURE FUNCTION GetMatrixColumns_l(this) RESULT(columns)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN) :: this
       INTEGER :: columns
     END FUNCTION GetMatrixColumns_l
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
     !! The current implementation has you go from matrix to triplet list,
     !! triplet list to transposed triplet list. The triplet list must then be
     !! sorted and then the return matrix is constructed.
     !! @param[in] this the matrix to be transposed.
     !! @param[out] matT the input matrix transposed.
     PURE SUBROUTINE TransposeMatrix_l(this, matT)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN)  :: this
       CLASS(Matrix_l), INTENT(INOUT) :: matT
     END SUBROUTINE TransposeMatrix_l
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
     !> Matrix B = alpha*Matrix A + Matrix B (AXPY).
     !! This will utilize the sparse vector addition routine.
     !! @param[in] matA Matrix A.
     !! @param[in,out] matB Matrix B.
     !! @param[in] alpha_in multiplier (optional, default=1.0)
     !! @param[in] threshold_in for flushing values to zero. (Optional, default=0).
     PURE SUBROUTINE IncrementMatrix_l(matA, matB, alpha_in, threshold_in)
       USE DataTypesModule, ONLY : NTREAL
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN)  :: matA
       CLASS(Matrix_l), INTENT(INOUT) :: matB
       REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
       REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
     END SUBROUTINE IncrementMatrix_l
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !> Pairwise Multiply two matrices.
     !! This will utilize the sparse vector pairwise routine.
     !! @param[in] matA Matrix A.
     !! @param[in] matB Matrix B.
     !! @param[in,out] matC = MatA mult MatB.
     PURE SUBROUTINE PairwiseMultiplyMatrix_l(matA, matB, matC)
       IMPORT :: Matrix_l
       IMPLICIT NONE
       !! Parameters
       CLASS(Matrix_l), INTENT(IN)  :: matA
       CLASS(Matrix_l), INTENT(IN) :: matB
       CLASS(Matrix_l), INTENT(INOUT) :: matC
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
       REAL(NTREAL), DIMENSION(this%columns), INTENT(OUT) :: norm_per_column
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixModule
