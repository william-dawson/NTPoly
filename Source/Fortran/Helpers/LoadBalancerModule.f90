!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for load balancing the matrix multiplication calculation.
MODULE LoadBalancerModule
  USE MatrixMemoryPoolPModule, ONLY : MatrixMemoryPool_p
  USE MatrixPSAlgebraModule, ONLY : DistributedGemm
  USE MatrixPSModule, ONLY : Matrix_ps, ConstructEmptyMatrixDS, &
       & DestructMatrixDS, FillDistributedPermutation
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
    TYPE(Matrix_ps), INTENT(IN) :: mat
    TYPE(Matrix_ps), INTENT(INOUT) :: mat_out
    TYPE(Permutation_t), INTENT(IN) :: permutation
    TYPE(MatrixMemoryPool_p), INTENT(INOUT),OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL ConstructEmptyDistributedSparseMatrix(PermuteRows, &
         & mat%actual_matrix_dimension, mat%process_grid)
    CALL ConstructEmptyDistributedSparseMatrix(PermuteColumns, &
         & mat%actual_matrix_dimension, mat%process_grid)
    CALL FillDistributedPermutation(PermuteRows, permutation%index_lookup, &
         & permuterows=.TRUE.)
    CALL FillDistributedPermutation(PermuteColumns, permutation%index_lookup, &
         & permuterows=.FALSE.)
    CALL ConstructEmptyDistributedSparseMatrix(Temp, &
         mat%actual_matrix_dimension, mat%process_grid)

    !! Permute Matrices.
    IF (PRESENT(memorypool_in)) THEN
       CALL DistributedGemm(PermuteRows, mat, Temp, &
            & memory_pool_in=memorypool_in)
       CALL DistributedGemm(Temp, PermuteColumns, mat_out, &
            & memory_pool_in=memorypool_in)
    ELSE
       CALL DistributedGemm(PermuteRows, mat, Temp)
       CALL DistributedGemm(Temp, PermuteColumns, mat_out)
    END IF

    CALL DestructDistributedSparseMatrix(PermuteRows)
    CALL DestructDistributedSparseMatrix(PermuteColumns)
    CALL DestructDistributedSparseMatrix(Temp)
  END SUBROUTINE PermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Undo a permutation applied to a matrix.
  !! @param[in] mat matrix to undo permutation of.
  !! @param[out] mat_out unpermuted matrix.
  !! @param[in] permutation to remove.
  !! @param[inout] memorypool_in memory pool to use. Optional.
  SUBROUTINE UndoPermuteMatrix(mat, mat_out, permutation, memorypool_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: mat
    TYPE(Matrix_ps), INTENT(INOUT) :: mat_out
    TYPE(Permutation_t), INTENT(IN) :: permutation
    TYPE(MatrixMemoryPool_p), INTENT(INOUT),OPTIONAL :: memorypool_in
    !! Local Variables
    TYPE(Matrix_ps) :: PermuteRows, PermuteColumns
    TYPE(Matrix_ps) :: Temp

    !! Build Permutation Matrices
    CALL ConstructEmptyDistributedSparseMatrix(PermuteRows, &
         & mat%actual_matrix_dimension, mat%process_grid)
    CALL ConstructEmptyDistributedSparseMatrix(PermuteColumns, &
         mat%actual_matrix_dimension, mat%process_grid)
    CALL FillDistributedPermutation(PermuteRows, permutation%index_lookup, &
         & permuterows=.TRUE.)
    CALL FillDistributedPermutation(PermuteColumns, permutation%index_lookup, &
         & permuterows=.FALSE.)
    CALL ConstructEmptyDistributedSparseMatrix(Temp, &
         & mat%actual_matrix_dimension, mat%process_grid)

    !! Permute Matrices.
    IF (PRESENT(memorypool_in)) THEN
       CALL DistributedGemm(PermuteColumns, mat, Temp, &
            & memory_pool_in=memorypool_in)
       CALL DistributedGemm(Temp, PermuteRows, mat_out, &
            & memory_pool_in=memorypool_in)
    ELSE
       CALL DistributedGemm(PermuteColumns, mat, Temp)
       CALL DistributedGemm(Temp, PermuteRows, mat_out)
    END IF

    !! Cleanup
    CALL DestructDistributedSparseMatrix(PermuteRows)
    CALL DestructDistributedSparseMatrix(PermuteColumns)
    CALL DestructDistributedSparseMatrix(Temp)
  END SUBROUTINE UndoPermuteMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE LoadBalancerModule
