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
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: DenseMatrix_t
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE DenseMatrix_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptyDenseMatrix
  PUBLIC :: ConstructDenseFromSparse
  PUBLIC :: ConstructSparseFromDense
  PUBLIC :: DestructDenseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MultiplyDense
  PUBLIC :: DenseEigenDecomposition
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  PURE SUBROUTINE ConstructEmptyDenseMatrix(this, rows, columns)
    !! Parameters
    TYPE(DenseMatrix_t), INTENT(INOUT) :: this
    INTEGER :: rows
    INTEGER :: columns

    CALL DestructDenseMatrix(this)
    ALLOCATE(this%data(rows,columns))

    this%rows = rows
    this%columns = columns

  END SUBROUTINE ConstructEmptyDenseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  !! @param[in] sparse_matrix a sparse matrix to convert.
  !! @param[inout] dense_matrix output. Must be preallocated.
  PURE SUBROUTINE ConstructDenseFromSparse(sparse_matrix, dense_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: sparse_matrix
    TYPE(DenseMatrix_t), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    INTEGER :: inner_counter, outer_counter
    INTEGER :: elements_per_inner
    INTEGER :: total_counter
    TYPE(Triplet_t) :: temporary

    IF (.NOT. ALLOCATED(dense_matrix%data)) THEN
       CALL ConstructEmptyDenseMatrix(dense_matrix,sparse_matrix%rows, &
            & sparse_matrix%columns)
    END IF

    !! Loop over elements.
    dense_matrix%data = 0
    total_counter = 1
    DO outer_counter = 1, sparse_matrix%columns
       elements_per_inner = sparse_matrix%outer_index(outer_counter+1) - &
            & sparse_matrix%outer_index(outer_counter)
       temporary%index_column = outer_counter
       DO inner_counter = 1, elements_per_inner
          temporary%index_row = sparse_matrix%inner_index(total_counter)
          temporary%point_value = sparse_matrix%values(total_counter)
          dense_matrix%data(temporary%index_row, temporary%index_column) = &
               & temporary%point_value
          total_counter = total_counter + 1
       END DO
    END DO
  END SUBROUTINE ConstructDenseFromSparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  !! @param[in] dense_matrix to convert.
  !! @param[out] sparse_matrix output matrix.
  !! @param[in] threshold_in value for pruning values to zero.
  PURE SUBROUTINE ConstructSparseFromDense(dense_matrix, sparse_matrix, &
       & threshold_in)
    !! Parameters
    TYPE(DenseMatrix_t), INTENT(IN) :: dense_matrix
    TYPE(SparseMatrix_t), INTENT(INOUT) :: sparse_matrix
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    INTEGER :: inner_counter, outer_counter
    TYPE(Triplet_t) :: temporary
    TYPE(TripletList_t) :: temporary_list
    INTEGER :: columns, rows
    INTEGER :: ind

    columns = dense_matrix%columns
    rows = dense_matrix%rows

    IF (PRESENT(threshold_in)) THEN
       CALL ConstructTripletList(temporary_list)
       DO outer_counter = 1, columns
          temporary%index_column = outer_counter
          DO inner_counter = 1, rows
             temporary%point_value = &
                  & dense_matrix%data(inner_counter,outer_counter)
             IF (ABS(temporary%point_value) .GT. threshold_in) THEN
                temporary%index_row = inner_counter
                CALL AppendToTripletList(temporary_list,temporary)
             END IF
          END DO
       END DO
    ELSE
       CALL ConstructTripletList(temporary_list, rows*columns)
       DO outer_counter = 1, columns
          temporary%index_column = outer_counter
          DO inner_counter = 1, rows
             temporary%point_value = &
                  & dense_matrix%data(inner_counter,outer_counter)
             temporary%index_row = inner_counter
             temporary_list%data(inner_counter+rows*(outer_counter-1)) = &
                  & temporary
          END DO
       END DO
    END IF

    CALL ConstructFromTripletList(sparse_matrix, temporary_list, rows, columns)
  END SUBROUTINE ConstructSparseFromDense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE DestructDenseMatrix(this)
    !! Parameters
    TYPE(DenseMatrix_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%data)) THEN
       DEALLOCATE(this%data)
    END IF
  END SUBROUTINE DestructDenseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  !! @param[in] MatA the first matrix.
  !! @param[in] MatB the second matrix.
  !! @param[inout] MatC = MatA*MatB.
  PURE SUBROUTINE MultiplyDense(MatA,MatB,MatC)
    !! Parameters
    TYPE(DenseMatrix_t), INTENT(IN) :: MatA
    TYPE(DenseMatrix_t), INTENT(IN) :: MatB
    TYPE(DenseMatrix_t), INTENT(INOUT) :: MatC

    IF (.NOT. ALLOCATED(MatC%data)) THEN
       CALL ConstructEmptyDenseMatrix(MatC,MatA%rows,MatB%columns)
    END IF

    MatC%data = MATMUL(MatA%data,MatB%data)
  END SUBROUTINE MultiplyDense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the eigenvectors of a dense matrix.
  !! Wraps a standard dense linear algebra routine.
  !! @param[in] MatA the matrix to decompose.
  !! @param[out] MatV the eigenvectors.
  SUBROUTINE DenseEigenDecomposition(MatA, MatV)
    !! Parameters
    TYPE(DenseMatrix_t), INTENT(IN) :: MatA
    TYPE(DenseMatrix_t), INTENT(INOUT) :: MatV
    !! Local variables
    CHARACTER, PARAMETER :: job = 'V', uplo = 'U'
    INTEGER :: N, LDA
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: W
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WORK
    DOUBLE PRECISION :: WORKTEMP
    INTEGER :: LWORK
    INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
    INTEGER :: IWORKTEMP
    INTEGER :: LIWORK
    INTEGER :: INFO

    IF (.NOT. ALLOCATED(MatV%data)) THEN
       CALL ConstructEmptyDenseMatrix(MatV,MatA%rows,MatA%columns)
    END IF
    MatV%data = MatA%data

    N = SIZE(MatA%data,DIM=1)
    LDA = N

    !! Allocations
    ALLOCATE(W(N))

    !! Determine the scratch space size
    LWORK = -1
    CALL DSYEVD(JOB, UPLO, N, MatA%data, LDA, W, WORKTEMP, LWORK, IWORKTEMP, &
         & LIWORK, INFO)
    N = LDA
    LWORK = INT(WORKTEMP)
    ALLOCATE(WORK(LWORK))
    LIWORK = INT(IWORKTEMP)
    ALLOCATE(IWORK(LIWORK))

    !! Run Lapack For Real
    CALL DSYEVD(JOB, UPLO, N, MatV%data, LDA, W, WORK, LWORK, IWORK, LIWORK, &
         & INFO)

    !! Cleanup
    DEALLOCATE(W)
    DEALLOCATE(Work)

  END SUBROUTINE DenseEigenDecomposition
END MODULE DenseMatrixModule
