!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for load balancing the matrix multiplication calculation.
MODULE LoadBalancerModule
  USE MatrixMemoryPoolPModule, ONLY : MatrixMemoryPool_p
  USE MatrixPSAlgebraModule, ONLY : MatrixMultiply
  USE MatrixPSModule, ONLY : Matrix_ps
  USE PermutationModule, ONLY : Permutation_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: PermuteMatrix
  PUBLIC :: UndoPermuteMatrix
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Apply a permutation to a matrix.
  !! @param[in] mat matrix to permute.
  !! @param[out] mat_out permuted matrix.
  !! @param[in] permutation to apply.
  !! @param[inout] memorypool_in memory pool to use. Optional.
  SUBROUTINE PermuteMatrix(mat, mat_out, permutation, memorypool_in)
    !! Parameters
    CLASS(Matrix_ps), INTENT(IN) :: mat
    CLASS(Matrix_ps), INTENT(INOUT) :: mat_out
    TYPE(Permutation_t), INTENT(IN) :: permutation
    CLASS(MatrixMemoryPool_p), INTENT(INOUT),OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL PermuteRows%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)
    CALL PermuteColumns%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)
    CALL PermuteRows%FillPermutation(permutation%index_lookup, &
         & permuterows=.TRUE.)
    CALL PermuteColumns%FillPermutation(permutation%index_lookup, &
         & permuterows=.FALSE.)
    CALL Temp%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)

    !! Permute Matrices.
    IF (PRESENT(memorypool_in)) THEN
       CALL MatrixMultiply(PermuteRows, mat, Temp, memory_pool_in=memorypool_in)
       CALL MatrixMultiply(Temp, PermuteColumns, mat_out, &
            & memory_pool_in=memorypool_in)
    ELSE
       CALL MatrixMultiply(PermuteRows, mat, Temp)
       CALL MatrixMultiply(Temp, PermuteColumns, mat_out)
    END IF

    !! Cleanup
    CALL PermuteRows%Destruct
    CALL PermuteColumns%Destruct
    CALL Temp%Destruct
  END SUBROUTINE PermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Undo a permutation applied to a matrix.
  !! @param[in] mat matrix to undo permutation of.
  !! @param[out] mat_out unpermuted matrix.
  !! @param[in] permutation to remove.
  !! @param[inout] memorypool_in memory pool to use. Optional.
  SUBROUTINE UndoPermuteMatrix(mat, mat_out, permutation, memorypool_in)
    !! Parameters
    CLASS(Matrix_ps), INTENT(IN) :: mat
    CLASS(Matrix_ps), INTENT(INOUT) :: mat_out
    TYPE(Permutation_t), INTENT(IN) :: permutation
    CLASS(MatrixMemoryPool_p), INTENT(INOUT),OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL PermuteRows%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)
    CALL PermuteColumns%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)
    CALL PermuteRows%FillPermutation(permutation%index_lookup, &
         & permuterows=.TRUE.)
    CALL PermuteColumns%FillPermutation(permutation%index_lookup, &
         & permuterows=.FALSE.)
    CALL Temp%InitEmpty(mat%actual_matrix_dimension, mat%process_grid)

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
    CALL PermuteRows%Destruct
    CALL PermuteColumns%Destruct
    CALL Temp%Destruct
  END SUBROUTINE UndoPermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoadBalancerModule
