!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for sending matrices between processes.
MODULE MatrixSendRecvModule
  USE DataTypesModule, ONLY : MPINTREAL
  USE MatrixDRModule, ONLY : Matrix_dr, ConstructEmptyMatrixD
  USE MatrixSRModule, ONLY : Matrix_sr, ConstructEmptyMatrixS
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a send call.
  TYPE, PUBLIC :: SendRecvHelper_t
     !> A request object for sending the sizes.
     INTEGER :: size_request
     !> A request object for sending outer indices.
     INTEGER :: outer_request
     !> A request object for sending inner indices.
     INTEGER :: inner_request
     !> A request object for sending data.
     INTEGER :: data_request
     !> Number of values to receive, rows, columns
     INTEGER, DIMENSION(3) :: matrix_data
  END TYPE SendRecvHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SendSparseMatrixSizes
  PUBLIC :: SendSparseMatrixData
  PUBLIC :: RecvSparseMatrixSizes
  PUBLIC :: RecvSparseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SendDenseMatrixSizes
  PUBLIC :: SendDenseMatrixData
  PUBLIC :: RecvDenseMatrixSizes
  PUBLIC :: RecvDenseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestSendRecvSizeRequest
  PUBLIC :: TestSendRecvOuterRequest
  PUBLIC :: TestSendRecvInnerRequest
  PUBLIC :: TestSendRecvDataRequest
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendSparseMatrixSizes(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    !! Local Data
    INTEGER :: ierr

    helper%matrix_data(1) = SIZE(inmat%values)
    helper%matrix_data(2) = inmat%rows
    helper%matrix_data(3) = inmat%columns

    CALL MPI_Isend(helper%matrix_data, 3, MPI_INTEGER, rank, send_tag, comm, &
         & helper%size_request, ierr)

  END SUBROUTINE SendSparseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvSparseMatrixSizes(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_sr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    !! Local Data
    INTEGER :: ierr

    CALL MPI_Irecv(helper%matrix_data, 3, MPI_INTEGER, rank, recv_tag, comm, &
         & helper%size_request, ierr)

  END SUBROUTINE RecvSparseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendSparseMatrixData(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Send the data
    CALL MPI_Isend(inmat%values, SIZE(inmat%values), MPINTREAL, &
         & rank, send_tag, comm, helper%data_request, ierr)
    CALL MPI_Isend(inmat%inner_index, SIZE(inmat%values), MPI_INT, &
         & rank, send_tag, comm, helper%inner_request, ierr)
    CALL MPI_Isend(inmat%outer_index, inmat%columns+1, MPI_INT, &
         & rank, send_tag, comm, helper%outer_request, ierr)

  END SUBROUTINE SendSparseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvSparseMatrixData(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_sr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    CALL ConstructEmptyMatrixS(outmat, helper%matrix_data(3), &
         & helper%matrix_data(2))
    ALLOCATE(outmat%values(helper%matrix_data(1)))
    ALLOCATE(outmat%inner_index(helper%matrix_data(1)))
    outmat%outer_index(1) = 0

    !! Receive the data
    CALL MPI_Irecv(outmat%values, helper%matrix_data(1), MPINTREAL, &
         & rank, recv_tag, comm, helper%data_request, ierr)
    CALL MPI_Irecv(outmat%inner_index, helper%matrix_data(1), MPI_INT, &
         & rank, recv_tag, comm, helper%inner_request, ierr)
    CALL MPI_Irecv(outmat%outer_index, helper%matrix_data(3)+1, MPI_INT, &
         & rank, recv_tag, comm, helper%outer_request, ierr)

  END SUBROUTINE RecvSparseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendDenseMatrixSizes(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_dr), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    !! Local Data
    INTEGER :: ierr

    helper%matrix_data(1) = inmat%rows*inmat%columns
    helper%matrix_data(2) = inmat%rows
    helper%matrix_data(3) = inmat%columns

    CALL MPI_Isend(helper%matrix_data, 3, MPI_INTEGER, rank, send_tag, comm, &
         & helper%size_request, ierr)

  END SUBROUTINE SendDenseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvDenseMatrixSizes(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_dr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    !! Local Data
    INTEGER :: ierr

    CALL MPI_Irecv(helper%matrix_data, 3, MPI_INTEGER, rank, recv_tag, comm, &
         & helper%size_request, ierr)

  END SUBROUTINE RecvDenseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendDenseMatrixData(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_dr), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Send the data
    CALL MPI_Isend(inmat%data, inmat%rows*inmat%columns, MPINTREAL, &
         & rank, send_tag, comm, helper%data_request, ierr)

  END SUBROUTINE SendDenseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvDenseMatrixData(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_dr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    CALL ConstructEmptyMatrixD(outmat, helper%matrix_data(3), &
         & helper%matrix_data(2))

    !! Receive the data
    CALL MPI_Irecv(outmat%data, helper%matrix_data(1), MPINTREAL, &
         & rank, recv_tag, comm, helper%data_request, ierr)

  END SUBROUTINE RecvDenseMatrixData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the send helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvSizeRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%size_request, request_completed, status, ierr)

  END FUNCTION TestSendRecvSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvOuterRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%outer_request, request_completed, status, ierr)

  END FUNCTION TestSendRecvOuterRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvInnerRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%inner_request, request_completed, status, ierr)

  END FUNCTION TestSendRecvInnerRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvDataRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%data_request, request_completed, status, ierr)

  END FUNCTION TestSendRecvDataRequest
END MODULE MatrixSendRecvModule
