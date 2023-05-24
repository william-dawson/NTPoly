!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains helper routines for converting an NTPoly matrix
!> to data structures used in other programs.
MODULE MatrixConversionModule
  USE DataTypesModule, ONLY : NTREAL
  USE PSMatrixModule, ONLY : Matrix_ps, ConvertMatrixToReal, CopyMatrix, &
       & DestructMatrix, MergeMatrixLocalBlocks, SplitMatrixToLocalBlocks
  USE PSMatrixAlgebraModule, ONLY : PairwiseMultiplyMatrix, ScaleMatrix, &
       & IncrementMatrix
  USE SMatrixModule, ONLY : Matrix_lsr, DestructMatrix
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
    TYPE(Matrix_lsr) :: local_mat

    !! First we need to make sure that the sparsity pattern is all 1s.
    IF (pattern%is_complex) THEN
       CALL ConvertMatrixToReal(pattern, pattern_1s)
    ELSE
       CALL CopyMatrix(pattern, pattern_1s)
    END IF
    CALL MergeMatrixLocalBlocks(pattern_1s, local_mat)
    local_mat%values = 1.0_NTREAL
    CALL SplitMatrixToLocalBlocks(pattern_1s, local_mat)

    !! Then all zeros
    CALL CopyMatrix(pattern_1s, pattern_0s)
    CALL ScaleMatrix(pattern_0s, 0.0_NTREAL)

    !! Here we add in the zero values that were missing from the original
    !! matrix. The secret here is that if you use a negative threshold, we
    !! never filter a value.
    CALL IncrementMatrix(pattern_0s, mat, threshold_in = -1.0_NTREAL)

    !! Next, we zero out values outside of the sparsity pattern.
    CALL CopyMatrix(mat, filtered)
    CALL PairwiseMultiplyMatrix(pattern_1s, filtered, mat)

    !! Cleanup
    CALL DestructMatrix(pattern_1s)
    CALL DestructMatrix(pattern_0s)
    CALL DestructMatrix(filtered)
    CALL DestructMatrix(local_mat)

  END SUBROUTINE SnapMatrixToSparsityPattern
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixConversionModule
