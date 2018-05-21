!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for broadcasting matrices across processes.
MODULE MatrixBroadcastModule
  USE DataTypesModule
  USE SparseMatrixModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a gather call.
  TYPE, PUBLIC :: BroadcastHelper_t
     !> A request object for broadcasting the sizes.
     INTEGER :: size_request
     !> A request object for broadcasting outer indices.
     INTEGER :: outer_request
     !> A request object for broadcasting inner indices.
     INTEGER :: inner_request
     !> A request object for broadcasting data.
     INTEGER :: data_request
     !> The rank of this process.
     INTEGER :: comm_rank
     !> Number of values to receive, rows, columns
     INTEGER, DIMENSION(3) :: matrix_data
  END TYPE BroadcastHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: BroadcastMatrixSizes
  PUBLIC :: BroadcastMatrixData
  PUBLIC :: TestBroadcastMatrixSize
  PUBLIC :: TestBroadcastMatrixOuter
  PUBLIC :: TestBroadcastMatrixInner
  PUBLIC :: TestBroadcastMatrixData
  PUBLIC :: BroadcastMatrix
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Broadcast The Size of a Matrix across the communicator.
  !! @param[inout] matrix to broadcast.
  !! @param[in] comm the communicator to broadcast along.
  !! @param[in] root the root process which holds the matrix.
  !! @param[inout] helper the broadcast matrix helper.
  SUBROUTINE BroadcastMatrixSizes(matrix, comm, root, helper)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, INTENT(INOUT) :: comm
    INTEGER, INTENT(IN) :: root
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: ierr

    CALL MPI_COMM_RANK(comm, helper%comm_rank, ierr)

    IF (helper%comm_rank .EQ. 0) THEN
       helper%matrix_data(1) = SIZE(matrix%values)
       helper%matrix_data(2) = matrix%rows
       helper%matrix_data(3) = matrix%columns
    ELSE
       CALL DestructSparseMatrix(matrix)
    END IF

    CALL MPI_Ibcast(helper%matrix_data, 3, MPI_INTEGER, root, comm, &
         & helper%size_request, ierr)

  END SUBROUTINE BroadcastMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Broadcast The Data of a Matrix across the communicator.
  !! @param[inout] matrix to broadcast.
  !! @param[in] comm the communicator to broadcast along.
  !! @param[in] root the root process which holds the matrix.
  !! @param[inout] helper the broadcast matrix helper.
  SUBROUTINE BroadcastMatrixData(matrix, comm, root, helper)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, INTENT(INOUT) :: comm
    INTEGER, INTENT(IN) :: root
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    IF (helper%comm_rank .NE. 0) THEN
       CALL ConstructEmptySparseMatrix(matrix, helper%matrix_data(3), &
            & helper%matrix_data(2))
       ALLOCATE(matrix%values(helper%matrix_data(1)))
       ALLOCATE(matrix%inner_index(helper%matrix_data(1)))
       matrix%outer_index(1) = 0
    END IF

    CALL MPI_Ibcast(matrix%values, helper%matrix_data(1), MPINTREAL, root, &
         & comm, helper%data_request, ierr)
    CALL MPI_Ibcast(matrix%inner_index, helper%matrix_data(1), MPI_INT, root, &
         & comm, helper%inner_request, ierr)
    CALL MPI_Ibcast(matrix%outer_index, helper%matrix_data(3)+1, MPI_INT, root,&
         & comm, helper%outer_request, ierr)

  END SUBROUTINE BroadcastMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the send helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestBroadcastMatrixSize(helper) RESULT(request_completed)
    !! Parameters
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%size_request, request_completed, status, ierr)

  END FUNCTION TestBroadcastMatrixSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestBroadcastMatrixOuter(helper) RESULT(request_completed)
    !! Parameters
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%outer_request, request_completed, status, ierr)

  END FUNCTION TestBroadcastMatrixOuter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestBroadcastMatrixInner(helper) RESULT(request_completed)
    !! Parameters
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%inner_request, request_completed, status, ierr)

  END FUNCTION TestBroadcastMatrixInner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestBroadcastMatrixData(helper) RESULT(request_completed)
    !! Parameters
    TYPE(BroadcastHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%data_request, request_completed, status, ierr)

  END FUNCTION TestBroadcastMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Broadcast a matrix from the root process across the communicator.
  !! @param[inout] matrix to broadcast.
  !! @param[in] comm the communicator to broadcast along.
  !! @param[in] root the root process which holds the matrix.
  SUBROUTINE BroadcastMatrix(matrix, comm, root)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matrix
    INTEGER, INTENT(INOUT) :: comm
    INTEGER, INTENT(IN) :: root
    !! Local Data
    INTEGER :: rank
    !! Matrix Buffer
    INTEGER :: matrows, matcolumns
    INTEGER :: nnz
    INTEGER :: ierr

    CALL MPI_COMM_RANK(comm, rank, ierr)

    !! Matrix Size
    IF (rank .EQ. root) THEN
       matrows = matrix%rows
       matcolumns = matrix%columns
       nnz = matrix%outer_index(matrix%columns+1)
    END IF
    CALL MPI_Bcast(matrows, 1, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(matcolumns, 1, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(nnz, 1, MPI_INT, root, comm, ierr)
    IF (rank .NE. root) THEN
       CALL ConstructZeroSparseMatrix(matrix, matrows, matcolumns)
       DEALLOCATE(matrix%inner_index)
       DEALLOCATE(matrix%values)
       ALLOCATE(matrix%inner_index(nnz))
       ALLOCATE(matrix%values(nnz))
    END IF

    !! Gather Matrix Data
    CALL MPI_Bcast(matrix%outer_index, matcolumns+1, MPI_INT, &
         & root, comm, ierr)
    CALL MPI_Bcast(matrix%inner_index, nnz, MPI_INT, root, comm, ierr)
    CALL MPI_Bcast(matrix%values, nnz, MPINTREAL, root, comm, ierr)
  END SUBROUTINE BroadcastMatrix
END MODULE MatrixBroadcastModule
