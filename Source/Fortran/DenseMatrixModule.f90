!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DenseMatrixModule
  USE DataTypesModule, ONLY : NTREAL
  USE SparseMatrixModule, ONLY : SparseMatrix_t, ConstructFromTripletList
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
       & AppendToTripletList
  USE TripletModule, ONLY : Triplet_t
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructDenseFromSparse
  PUBLIC :: ConstructSparseFromDense
  PUBLIC :: MultiplyDense
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  !! @param[in] sparse_matrix a sparse matrix to convert.
  !! @param[inout] dense_matrix output. Must be preallocated.
  PURE SUBROUTINE ConstructDenseFromSparse(sparse_matrix, dense_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: sparse_matrix
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: dense_matrix
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
          dense_matrix(temporary%index_row, temporary%index_column) = &
               & temporary%point_value
          total_counter = total_counter + 1
       END DO
    END DO
  END SUBROUTINE ConstructDenseFromSparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  !! @param[in] dense_matrix to convert.
  !! @param[out] sparse_matrix output matrix.
  !! @param[in] threshold value for pruning values to zero.
  PURE SUBROUTINE ConstructSparseFromDense(dense_matrix,sparse_matrix,threshold)
    !! Parameters
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: dense_matrix
    TYPE(SparseMatrix_t), INTENT(INOUT) :: sparse_matrix
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    INTEGER :: inner_counter, outer_counter
    TYPE(Triplet_t) :: temporary
    TYPE(TripletList_t) :: temporary_list
    INTEGER :: columns, rows

    columns = SIZE(dense_matrix,DIM=2)
    rows = SIZE(dense_matrix,DIM=1)

    CALL ConstructTripletList(temporary_list)
    DO outer_counter = 1, columns
       temporary%index_column = outer_counter
       DO inner_counter = 1, rows
          temporary%point_value = dense_matrix(inner_counter,outer_counter)
          IF (ABS(temporary%point_value) .GT. threshold) THEN
             temporary%index_row = inner_counter
             CALL AppendToTripletList(temporary_list,temporary)
          END IF
       END DO
    END DO

    CALL ConstructFromTripletList(sparse_matrix, temporary_list, rows, columns)
  END SUBROUTINE ConstructSparseFromDense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  !! @param[in] MatA the first matrix.
  !! @param[in] MatB the second matrix.
  !! @param[inout] MatC = MatA*MatB. MatC must be preallocated.
  PURE SUBROUTINE MultiplyDense(MatA,MatB,MatC)
    !! Parameters
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: MatA
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: MatB
    REAL(NTREAL), DIMENSION(:,:), INTENT(INOUT) :: MatC

    MatC = MATMUL(MatA,MatB)
  END SUBROUTINE MultiplyDense
END MODULE DenseMatrixModule
