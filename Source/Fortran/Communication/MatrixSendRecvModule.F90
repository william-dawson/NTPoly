!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for sending matrices between processes.
MODULE MatrixSendRecvModule
  USE DataTypesModule, ONLY : MPINTREAL, MPINTCOMPLEX
  USE MatrixDModule, ONLY : Matrix_ldr, Matrix_ldc, ConstructEmptyMatrix
  USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc, ConstructEmptyMatrix
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
     !> If the matrix being sent is dense
     LOGICAL :: is_dense
     !> Keep track of what data has been received
     INTEGER :: recv_stage
  END TYPE SendRecvHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SendMatrixSizes
  PUBLIC :: SendMatrixData
  PUBLIC :: RecvMatrixSizes
  PUBLIC :: RecvMatrixData
  PUBLIC :: TestSendRecvSizes
  PUBLIC :: TestSendRecvData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE SendMatrixSizes
     MODULE PROCEDURE SendMatrixSizes_lsr
     MODULE PROCEDURE SendMatrixSizes_lsc
     MODULE PROCEDURE SendMatrixSizes_ldr
     MODULE PROCEDURE SendMatrixSizes_ldc
  END INTERFACE
  INTERFACE SendMatrixData
     MODULE PROCEDURE SendMatrixData_lsr
     MODULE PROCEDURE SendMatrixData_lsc
     MODULE PROCEDURE SendMatrixData_ldr
     MODULE PROCEDURE SendMatrixData_ldc
  END INTERFACE
  INTERFACE RecvMatrixSizes
     MODULE PROCEDURE RecvMatrixSizes_lsr
     MODULE PROCEDURE RecvMatrixSizes_lsc
     MODULE PROCEDURE RecvMatrixSizes_ldr
     MODULE PROCEDURE RecvMatrixSizes_ldc
  END INTERFACE
  INTERFACE RecvMatrixData
     MODULE PROCEDURE RecvMatrixData_lsr
     MODULE PROCEDURE RecvMatrixData_lsc
     MODULE PROCEDURE RecvMatrixData_ldr
     MODULE PROCEDURE RecvMatrixData_ldc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixSizes_lsr(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixSizes.f90"

    helper%is_dense = .FALSE.

  END SUBROUTINE SendMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixSizes_lsc(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixSizes.f90"

    helper%is_dense = .FALSE.

  END SUBROUTINE SendMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixSizes_ldr(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixSizes.f90"

    helper%is_dense = .TRUE.

  END SUBROUTINE SendMatrixSizes_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send size information about matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixSizes_ldc(inmat, rank, comm, helper, send_tag)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixSizes.f90"

    helper%is_dense = .TRUE.

  END SUBROUTINE SendMatrixSizes_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixSizes_lsr(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: outmat

    INCLUDE "includes/RecvMatrixSizes.f90"

    helper%is_dense = .FALSE.

  END SUBROUTINE RecvMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixSizes_lsc(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: outmat

    INCLUDE "includes/RecvMatrixSizes.f90"

    helper%is_dense = .FALSE.

  END SUBROUTINE RecvMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixSizes_ldr(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    !! Local Data
    INTEGER :: ierr

    CALL MPI_Irecv(helper%matrix_data, 3, MPI_INTEGER, rank, recv_tag, comm, &
         & helper%size_request, ierr)

    helper%is_dense = .TRUE.

  END SUBROUTINE RecvMatrixSizes_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive size information about matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixSizes_ldc(outmat, rank, comm, helper, recv_tag)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    INTEGER, INTENT(INOUT) :: comm
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    !! Local Data
    INTEGER :: ierr

    CALL MPI_Irecv(helper%matrix_data, 3, MPI_INTEGER, rank, recv_tag, comm, &
         & helper%size_request, ierr)

    helper%is_dense = .TRUE.

  END SUBROUTINE RecvMatrixSizes_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixData_lsr(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixData.f90"
    CALL MPI_Isend(inmat%values, SIZE(inmat%values), MPINTREAL, &
         & rank, send_tag, comm, helper%data_request, ierr)

  END SUBROUTINE SendMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixData_lsc(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN) :: inmat

    INCLUDE "includes/SendMatrixData.f90"
    CALL MPI_Isend(inmat%values, SIZE(inmat%values), MPINTCOMPLEX, &
         & rank, send_tag, comm, helper%data_request, ierr)

  END SUBROUTINE SendMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixData_ldr(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Send the data
    CALL MPI_Isend(inmat%data, inmat%rows*inmat%columns, MPINTREAL, &
         & rank, send_tag, comm, helper%data_request, ierr)

  END SUBROUTINE SendMatrixData_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Send data contained in matrices.
  !! @param[in] inmat the matrix to send.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] send_tag a tag to identify this message.
  SUBROUTINE SendMatrixData_ldc(inmat,rank,comm,helper,send_tag)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(IN) :: inmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: send_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Send the data
    CALL MPI_Isend(inmat%data, inmat%rows*inmat%columns, MPINTCOMPLEX, &
         & rank, send_tag, comm, helper%data_request, ierr)

  END SUBROUTINE SendMatrixData_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixData_lsr(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(INOUT) :: outmat

    INCLUDE "includes/RecvMatrixData.f90"
    CALL MPI_Irecv(outmat%values, helper%matrix_data(1), MPINTREAL, &
         & rank, recv_tag, comm, helper%data_request, ierr)

  END SUBROUTINE RecvMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixData_lsc(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(INOUT) :: outmat

    INCLUDE "includes/RecvMatrixData.f90"
    CALL MPI_Irecv(outmat%values, helper%matrix_data(1), MPINTCOMPLEX, &
         & rank, recv_tag, comm, helper%data_request, ierr)

  END SUBROUTINE RecvMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixData_ldr(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    outmat = Matrix_ldr(helper%matrix_data(3), helper%matrix_data(2))

    !! Receive the data
    CALL MPI_Irecv(outmat%data, helper%matrix_data(1), MPINTREAL, &
         & rank, recv_tag, comm, helper%data_request, ierr)

  END SUBROUTINE RecvMatrixData_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Receive data contained in matrices.
  !! @param[inout] outmat the matrix to receive.
  !! @param[in] rank the rank of the current process.
  !! @param[inout] comm the communicator to perform sends along.
  !! @param[inout] helper the send-receive helper data structure.
  !! @param[in] recv_tag a tag to identify this message.
  SUBROUTINE RecvMatrixData_ldc(outmat,rank,comm,helper,recv_tag)
    !! Parameters
    TYPE(Matrix_ldc), INTENT(INOUT) :: outmat
    INTEGER, INTENT(IN) :: rank
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    INTEGER, INTENT(IN) :: recv_tag
    INTEGER :: comm
    !! Local Data
    INTEGER :: ierr

    !! Build Storage
    outmat = Matrix_ldc(helper%matrix_data(3), helper%matrix_data(2))

    !! Receive the data
    CALL MPI_Irecv(outmat%data, helper%matrix_data(1), MPINTCOMPLEX, &
         & rank, recv_tag, comm, helper%data_request, ierr)

  END SUBROUTINE RecvMatrixData_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the send helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvSizes(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    CALL MPI_Test(helper%size_request, request_completed, status, ierr)

  END FUNCTION TestSendRecvSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the matrix data is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestSendRecvData(helper) RESULT(request_completed)
    !! Parameters
    TYPE(SendRecvHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed
    !! Local Data
    INTEGER :: ierr
    INTEGER :: status(MPI_STATUS_SIZE)

    IF (helper%is_dense) THEN
       IF (helper%recv_stage .EQ. 0) THEN
          CALL MPI_Test(helper%data_request, request_completed, status, ierr)
          IF (request_completed) THEN
            helper%recv_stage = 1
          END IF
       ELSE
          request_completed = .TRUE.
       END IF
    ELSE
       IF (helper%recv_stage .EQ. 0) THEN
          CALL MPI_Test(helper%outer_request, request_completed, status, ierr)
          IF (request_completed) THEN
            helper%recv_stage = 1
          END IF
       ELSE IF (helper%recv_stage .EQ. 1) THEN
          CALL MPI_Test(helper%inner_request, request_completed, status, ierr)
          IF (request_completed) THEN
            helper%recv_stage = 2
          END IF
       ELSE IF (helper%recv_stage .EQ. 2) THEN
          CALL MPI_Test(helper%data_request, request_completed, status, ierr)
          IF (request_completed) THEN
            helper%recv_stage = 3
          END IF
       ELSE
          request_completed = .TRUE.
       END IF
    END IF

  END FUNCTION TestSendRecvData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixSendRecvModule
