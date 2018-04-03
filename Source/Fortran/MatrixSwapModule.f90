!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for swapping matrices between processes.
MODULE MatrixSwapModule
  USE DataTypesModule
  USE SparseMatrixModule
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a send call.
  TYPE, PUBLIC :: SwapHelper_t
     !> A request object for sending the sizes.
     INTEGER, DIMENSION(2) :: size_request
     !> A request object for sending outer indices.
     INTEGER, DIMENSION(2) :: outer_request
     !> A request object for sending inner indices.
     INTEGER, DIMENSION(2) :: inner_request
     !> A request object for sending data.
     INTEGER, DIMENSION(2) :: data_request
     !> Number of values to receive
     INTEGER :: num_values
  END TYPE SwapHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SwapMatrixSizes
  PUBLIC :: SwapMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestSwapSizeRequest
  PUBLIC :: TestSwapOuterRequest
  PUBLIC :: TestSwapInnerRequest
  PUBLIC :: TestSwapDataRequest
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Swap size information about matrices.
  SUBROUTINE SwapMatrixSizes(inmat, partner_rank, comm, helper, tag)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: partner_rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: tag
    !! Local Data
    INTEGER :: ierr

    CALL MPI_Isend(inmat%values, 1, MPI_INTEGER, partner_rank, tag, comm, &
         & helper%size_request(1), ierr)
    CALL MPI_Irecv(helper%num_values, 1, MPI_INTEGER, partner_rank, tag, comm, &
         & helper%size_request(2), ierr)
  END SUBROUTINE SwapMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Swap data contained in matrices.
  SUBROUTINE SwapMatrixData(inmat, outmat, partner_rank, comm, helper, tag)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: inmat
    TYPE(SparseMatrix_t), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: partner_rank
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    CALL ConstructEmptySparseMatrix(outmat, inmat%columns, inmat%rows)
    ALLOCATE(outmat%values(helper%num_values))
    ALLOCATE(outmat%inner_index(helper%num_values))
    outmat%outer_index(1) = 0

    !! Send/Recv the data
    CALL MPI_Isend(inmat%values, SIZE(inmat%values), MPINTREAL, &
         & partner_rank, tag, comm, helper%data_request(1), ierr)
    CALL MPI_Isend(inmat%inner_index, SIZE(inmat%values), MPI_INT, &
         & partner_rank, tag, comm, helper%inner_request(1), ierr)
    CALL MPI_Isend(inmat%outer_index, inmat%columns+1, MPI_INT, &
         & partner_rank, tag, comm, helper%outer_request(1), ierr)

    CALL MPI_Irecv(outmat%values, helper%num_values, MPINTREAL, &
         & partner_rank, tag, comm, helper%data_request(2), ierr)
    CALL MPI_Irecv(outmat%inner_index, helper%num_values, MPI_INT, &
         & partner_rank, tag, comm, helper%inner_request(2), ierr)
    CALL MPI_Irecv(outmat%outer_index, inmat%columns+1, MPI_INT, &
         & partner_rank, tag, comm, helper%outer_request(2), ierr)

  END SUBROUTINE SwapMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the send helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSwapSizeRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    LOGICAL :: send_size, recv_size
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%size_request(1), send_size, status, ierr)
    CALL MPI_Test(helper%size_request(2), recv_size, status, ierr)

    request_completed = send_size .AND. recv_size
  END FUNCTION TestSwapSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSwapOuterRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    LOGICAL :: send_size, recv_size
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%outer_request(1), send_size, status, ierr)
    CALL MPI_Test(helper%outer_request(2), recv_size, status, ierr)

    request_completed = send_size .AND. recv_size
  END FUNCTION TestSwapOuterRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSwapInnerRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    LOGICAL :: send_size, recv_size
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%inner_request(1), send_size, status, ierr)
    CALL MPI_Test(helper%inner_request(2), recv_size, status, ierr)

    request_completed = send_size .AND. recv_size
  END FUNCTION TestSwapInnerRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSwapDataRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SwapHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    LOGICAL :: send_size, recv_size
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%data_request(1), send_size, status, ierr)
    CALL MPI_Test(helper%data_request(2), recv_size, status, ierr)

    request_completed = send_size .AND. recv_size
  END FUNCTION TestSwapDataRequest
END MODULE MatrixSwapModule
