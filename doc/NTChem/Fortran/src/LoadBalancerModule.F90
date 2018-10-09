!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for load balancing the matrix multiplication calculation.
MODULE LoadBalancerModule
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & DestructMatrix, FillMatrixPermutation
  USE PermutationModule, ONLY : Permutation_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PermuteMatrix
  PUBLIC :: UndoPermuteMatrix
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Apply a permutation to a matrix.
  SUBROUTINE PermuteMatrix(mat, mat_out, permutation, memorypool_in)
    !> The matrix to permute.
    TYPE(Matrix_ps), INTENT(IN) :: mat
    !> The permuted matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: mat_out
    !> The permutation to apply.
    TYPE(Permutation_t), INTENT(IN) :: permutation
    !> Memory pool to use
    TYPE(MatrixMemoryPool_p), INTENT(INOUT), OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL ConstructEmptyMatrix(PermuteRows, mat)
    CALL ConstructEmptyMatrix(PermuteColumns, mat)
    CALL FillMatrixPermutation(PermuteRows, permutation%index_lookup, &
         & permute_rows_in=.TRUE.)
    CALL FillMatrixPermutation(PermuteColumns, permutation%index_lookup, &
         & permute_rows_in=.FALSE.)
    CALL ConstructEmptyMatrix(Temp, mat)

    !! Permute Matrices.
    IF (PRESENT(memorypool_in)) THEN
       CALL MatrixMultiply(PermuteRows, mat, Temp, &
            & memory_pool_in=memorypool_in)
       CALL MatrixMultiply(Temp, PermuteColumns, mat_out, &
            & memory_pool_in=memorypool_in)
    ELSE
       CALL MatrixMultiply(PermuteRows, mat, Temp)
       CALL MatrixMultiply(Temp, PermuteColumns, mat_out)
    END IF

    CALL DestructMatrix(PermuteRows)
    CALL DestructMatrix(PermuteColumns)
    CALL DestructMatrix(Temp)
  END SUBROUTINE PermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Undo a permutation applied to a matrix.
  SUBROUTINE UndoPermuteMatrix(mat, mat_out, permutation, memorypool_in)
    !> Matrix to undo permutation of.
    TYPE(Matrix_ps), INTENT(IN) :: mat
    !> Unpermuted matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: mat_out
    !> Permutation to remove.
    TYPE(Permutation_t), INTENT(IN) :: permutation
    !> Memory pool to use.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT),OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL ConstructEmptyMatrix(PermuteRows, mat)
    CALL ConstructEmptyMatrix(PermuteColumns, mat)
    CALL FillMatrixPermutation(PermuteRows, permutation%index_lookup, &
         & permute_rows_in=.TRUE.)
    CALL FillMatrixPermutation(PermuteColumns, permutation%index_lookup, &
         & permute_rows_in=.FALSE.)
    CALL ConstructEmptyMatrix(Temp, mat)

    !! Permute Matrices.
    IF (PRESENT(memorypool_in)) THEN
       CALL MatrixMultiply(PermuteColumns, mat, Temp, &
            & memory_pool_in=memorypool_in)
       CALL MatrixMultiply(Temp, PermuteRows, mat_out, &
            & memory_pool_in=memorypool_in)
    ELSE
       CALL MatrixMultiply(PermuteColumns, mat, Temp)
       CALL MatrixMultiply(Temp, PermuteRows, mat_out)
    END IF

    !! Cleanup
    CALL DestructMatrix(PermuteRows)
    CALL DestructMatrix(PermuteColumns)
    CALL DestructMatrix(Temp)
  END SUBROUTINE UndoPermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoadBalancerModule
