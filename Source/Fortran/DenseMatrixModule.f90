!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DenseMatrixModule
  USE DataTypesModule
  USE SparseMatrixModule
  USE TripletListModule
  USE TripletModule
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ConstructDenseFromSparse
  PUBLIC :: ConstructSparseFromDense
  PUBLIC :: MultiplyDense
  PUBLIC :: DenseEigenDecomposition
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  !! @param[in] threshold for pruning values to zero.
  SUBROUTINE DenseEigenDecomposition(MatA, MatV, threshold)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(In) :: MatA
    TYPE(SparseMatrix_t), INTENT(INOUT) :: MatV
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Matrices
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DMatA
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DMatV
    !! Local variables
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: A
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
    INTEGER :: LWORK
    DOUBLE PRECISION :: TEMP
    INTEGER :: INFO
    !! Temporary variables
    INTEGER :: II, JJ, counter

    !! Convert Input To Dense
    ALLOCATE(DMatA(MatA%rows, MatA%columns))
    DMatA = 0
    ALLOCATE(DMatV(MatA%rows, MatA%columns))
    DMatV = 0
    CALL ConstructDenseFromSparse(MatA,DMatA)

    N = SIZE(DMatA,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(A(N*N))
    ALLOCATE(W(N))

    !! Store as an upper triangular matrix
    A = 0
    counter = 1
    DO II = 1, N
       DO JJ = 1, N
          A(counter) = DMatA(II,JJ)
          counter = counter + 1
       END DO
    END DO

    !! Determine the scratch space size
    LWORK = -1
    CALL DSYEV(JOB, UPLO, N, A, LDA, W, TEMP, LWORK, INFO)
    N = LDA
    LWORK = INT(TEMP)
    ALLOCATE(WORK(LWORK))

    !! Run Lapack For Real
    CALL DSYEV(JOB, UPLO, N, A, LDA, W, WORK, LWORK, INFO)

    !! Unpack
    counter = 1
    DO II = 1, N
       DO JJ = 1, N
          DMatV(JJ,II) = A(counter)
          counter = counter + 1
       END DO
    END DO

    !! Convert Output To Sparse
    CALL ConstructSparseFromDense(DMatV,MatV,threshold)

    !! Cleanup
    DEALLOCATE(DMatA)
    DEALLOCATE(DMatV)
    DEALLOCATE(A)
    DEALLOCATE(W)
    DEALLOCATE(Work)

  END SUBROUTINE DenseEigenDecomposition
END MODULE DenseMatrixModule
