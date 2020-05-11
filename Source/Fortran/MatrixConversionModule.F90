!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains helper routines for converting an NTPoly matrix
!> to data structures used in other programs.
MODULE MatrixConversionModule
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixMapsModule, ONLY : MapMatrix_psr
  USE PSMatrixModule, ONLY : Matrix_ps, ConvertMatrixToReal, CopyMatrix, &
       & DestructMatrix
  USE PSMatrixAlgebraModule, ONLY : PairwiseMultiplyMatrix, ScaleMatrix, &
       & IncrementMatrix
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SnapMatrixToSparsityPattern
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Some codes use a fixed sparsity pattern for a matrix instead of filtering
  !> small values. Using this routine, the matrix is filled to have the same
  !> pattern as the second matrix argument. Zeros of the sparsity pattern are
  !> left in, whereas values outside the sparsity are removed. This can 
  !> faciliate conversion between formats.
  SUBROUTINE SnapMatrixToSparsityPattern(mat, pattern)
    !> The matrix to modify.
    TYPE(Matrix_ps), INTENT(INOUT) :: mat
    !> The matrix which defines the sparsity pattern.
    TYPE(Matrix_ps), INTENT(IN) :: pattern
    !! Local Variables
    TYPE(Matrix_ps) :: filtered
    TYPE(Matrix_ps) :: pattern_1s
    TYPE(Matrix_ps) :: pattern_0s
    TYPE(Matrix_ps) :: pattern_real
    INTEGER :: II

    !! First we need to make sure that the sparsity pattern is all 1s.
    IF (pattern%is_complex) THEN
       CALL ConvertMatrixToReal(pattern, pattern_real)
    ELSE
       CALL CopyMatrix(pattern, pattern_real)
    END IF
    CALL MapMatrix_psr(pattern_real, pattern_1s, SetMatrixToOne)

    !! Then all zeros
    CALL CopyMatrix(pattern_1s, pattern_0s)
    CALL ScaleMatrix(pattern_0s, 0.0_NTREAL)

    !! Here we add in the zero values that were missing from the original
    !! matrix. The secret here is that if you use a negative threshold, we
    !! never filter a value.
    CALL IncrementMatrix(pattern_0s, mat, threshold_in=-1.0_NTREAL)

    !! Next, we zero out values outside of the sparsity pattern.
    CALL CopyMatrix(mat, filtered)
    CALL PairwiseMultiplyMatrix(pattern_1s, filtered, mat)

    !! Cleanup
    CALL DestructMatrix(pattern_1s)
    CALL DestructMatrix(pattern_0s)
    CALL DestructMatrix(pattern_real)
    CALL DestructMatrix(filtered)

  END SUBROUTINE SnapMatrixToSparsityPattern
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION SetMatrixToOne(row, column, val) RESULT(valid)
    INTEGER, INTENT(INOUT), OPTIONAL :: row
    INTEGER, INTENT(INOUT), OPTIONAL :: column
    REAL(NTREAL), INTENT(INOUT), OPTIONAL :: val
    LOGICAL :: valid

    val = 1.0_NTREAL
    valid = .TRUE.
  END FUNCTION SetMatrixToOne
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE