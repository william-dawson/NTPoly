!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This (currently under construction) module allows one to convert
!! a sparse matrix to a dense matrix. It also supports dense the dense
!! versions of core matrix routines. This module will be used in situations
!! where matrices become too dense for good sparse matrix performance.
MODULE DenseMatrixModule
  USE DataTypesModule
  USE SparseMatrixModule
  USE TripletModule
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ConstructDenseFromSparse
  PUBLIC :: MultiplyDense
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  !! @param[in] sparse_matrix a sparse matrix to convert.
  !! @param[inout] dense_matrix a preallocated
  PURE SUBROUTINE ConstructDenseFromSparse(sparse_matrix,dense_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: sparse_matrix
    REAL(NTREAL), DIMENSION(:,:), INTENT(inout) :: dense_matrix
    !! Helper Variables
    INTEGER :: inner_counter, outer_counter
    INTEGER :: elements_per_inner
    INTEGER :: total_counter
    TYPE(Triplet_t) :: temporary

    !! Loop over elements.
    dense_matrix = 0
    total_counter = 1
    DO outer_counter = 1, sparse_matrix%columns
       elements_per_inner = sparse_matrix%outer_index(outer_counter+1) - &
            & sparse_matrix%outer_index(outer_counter)
       temporary%index_column = outer_counter
       DO inner_counter = 1, elements_per_inner
          temporary%index_row = sparse_matrix%inner_index(total_counter)
          temporary%point_value = sparse_matrix%values(total_counter)
          dense_matrix(temporary%index_column,temporary%index_row) = &
               & temporary%point_value
          total_counter = total_counter + 1
       END DO
    END DO
  END SUBROUTINE ConstructDenseFromSparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  !! @param[in] MatA the first matrix.
  !! @param[in] MatB the second matrix.
  !! @param[inout] MatC = MatA*MatB
  PURE SUBROUTINE MultiplyDense(MatA,MatB,MatC)
    !! Parameters
    REAL(NTREAL), DIMENSION(:,:), INTENT(in) :: MatA
    REAL(NTREAL), DIMENSION(:,:), INTENT(in) :: MatB
    REAL(NTREAL), DIMENSION(:,:), INTENT(inout) :: MatC

    MatC = MATMUL(MatA,MatB)
  END SUBROUTINE MultiplyDense
END MODULE DenseMatrixModule
