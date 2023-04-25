!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for reducing matrices across processes.
MODULE MatrixReduceModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, MPINTCOMPLEX, MPINTINTEGER
  USE SMatrixAlgebraModule, ONLY : IncrementMatrix
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, ConstructEmptyMatrix, &
       & DestructMatrix, CopyMatrix
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a reduce call.
  TYPE, PUBLIC :: ReduceHelper_t
     !> Number of processors involved in this gather.
     INTEGER :: comm_size
     !> A request object for gathering outer indices.
     INTEGER :: outer_request
     !> A request object for gathering inner indices.
     INTEGER :: inner_request
     !> A request object for gathering data.
     INTEGER :: data_request
     !> The error code after an MPI call.
     INTEGER :: error_code
     !> Number of values to gather from each process.
     INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_process
     !> The displacements for where those gathered values should go.
     INTEGER, DIMENSION(:), ALLOCATABLE :: displacement
#ifdef NOIALLGATHER
     !> For mpi backup, a list of request objets for outer indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_send_request_list
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_recv_request_list
     !> For mpi backup, a list of request objects for inner indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_send_request_list
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_recv_request_list
     !> For mpi backup, a list of request object for data.
     INTEGER, DIMENSION(:), ALLOCATABLE :: data_send_request_list
     INTEGER, DIMENSION(:), ALLOCATABLE :: data_recv_request_list
#endif
  END TYPE ReduceHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndComposeMatrixSizes
  PUBLIC :: ReduceAndComposeMatrixData
  PUBLIC :: ReduceAndComposeMatrixCleanup
  PUBLIC :: ReduceAndComposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndSumMatrixSizes
  PUBLIC :: ReduceAndSumMatrixData
  PUBLIC :: ReduceAndSumMatrixCleanup
  PUBLIC :: ReduceAndSumMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestReduceSizeRequest
  PUBLIC :: TestReduceInnerRequest
  PUBLIC :: TestReduceDataRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ReduceAndComposeMatrixSizes
     MODULE PROCEDURE ReduceAndComposeMatrixSizes_lsr
     MODULE PROCEDURE ReduceAndComposeMatrixSizes_lsc
  END INTERFACE ReduceAndComposeMatrixSizes
  INTERFACE ReduceAndComposeMatrixData
     MODULE PROCEDURE ReduceAndComposeMatrixData_lsr
     MODULE PROCEDURE ReduceAndComposeMatrixData_lsc
  END INTERFACE ReduceAndComposeMatrixData
  INTERFACE ReduceAndComposeMatrixCleanup
     MODULE PROCEDURE ReduceAndComposeMatrixCleanup_lsr
     MODULE PROCEDURE ReduceAndComposeMatrixCleanup_lsc
  END INTERFACE ReduceAndComposeMatrixCleanup
  INTERFACE ReduceAndComposeMatrix
     MODULE PROCEDURE ReduceAndComposeMatrix_lsr
     MODULE PROCEDURE ReduceAndComposeMatrix_lsc
  END INTERFACE ReduceAndComposeMatrix
  INTERFACE ReduceAndSumMatrixSizes
     MODULE PROCEDURE ReduceAndSumMatrixSizes_lsr
     MODULE PROCEDURE ReduceAndSumMatrixSizes_lsc
  END INTERFACE ReduceAndSumMatrixSizes
  INTERFACE ReduceAndSumMatrixData
     MODULE PROCEDURE ReduceAndSumMatrixData_lsr
     MODULE PROCEDURE ReduceAndSumMatrixData_lsc
  END INTERFACE ReduceAndSumMatrixData
  INTERFACE ReduceAndSumMatrixCleanup
     MODULE PROCEDURE ReduceAndSumMatrixCleanup_lsr
     MODULE PROCEDURE ReduceAndSumMatrixCleanup_lsc
  END INTERFACE ReduceAndSumMatrixCleanup
  INTERFACE ReduceAndSumMatrix
     MODULE PROCEDURE ReduceAndSumMatrix_lsr
     MODULE PROCEDURE ReduceAndSumMatrix_lsc
  END INTERFACE ReduceAndSumMatrix
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceAndComposeMatrixSizes_lsr(matrix, comm, gathered_matrix, &
    & helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The  helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndComposeMatrixSizes_sendrecv.f90"
#else
#include "comm_includes/ReduceAndComposeMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceAndComposeMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceAndComposeMatrixSizes_lsc(matrix, comm, gathered_matrix, &
    & helper)
    !! The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !! The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !! The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndComposeMatrixSizes_sendrecv.f90"
#else
#include "comm_includes/ReduceAndComposeMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceAndComposeMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align one matrix
  !> next to another.
  SUBROUTINE ReduceAndComposeMatrixData_lsr(matrix, comm, gathered_matrix, &
    & helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndComposeMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTREAL, &
            & II - 1, 4, comm, helper%data_send_request_list(II), ierr)
       istart = helper%displacement(II) + 1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, MPINTREAL, &
            & II - 1, 4, comm, helper%data_recv_request_list(II), ierr)
    END DO
#else
#include "comm_includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, comm, helper%data_request, ierr)
#endif
  END SUBROUTINE ReduceAndComposeMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align one matrix
  !> next to another.
  SUBROUTINE ReduceAndComposeMatrixData_lsc(matrix, comm, gathered_matrix, &
    & helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndComposeMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTCOMPLEX, &
            & II - 1, 4, comm, helper%data_send_request_list(II), ierr)
       istart = helper%displacement(II) + 1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, &
            & MPINTCOMPLEX, II - 1, 4, comm, &
            & helper%data_recv_request_list(II), ierr)
    END DO
#else
#include "comm_includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, comm, helper%data_request, ierr)
#endif
  END SUBROUTINE ReduceAndComposeMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third function to call, finishes setting up the matrices.
  PURE SUBROUTINE ReduceAndComposeMatrixCleanup_lsr(matrix, gathered_matrix, &
       & helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT) :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper

#include "comm_includes/ReduceAndComposeMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%outer_send_request_list)) THEN
       DEALLOCATE(helper%outer_send_request_list)
    END IF
    IF (ALLOCATED(helper%outer_recv_request_list)) THEN
       DEALLOCATE(helper%outer_recv_request_list)
    END IF
    IF (ALLOCATED(helper%inner_send_request_list)) THEN
       DEALLOCATE(helper%inner_send_request_list)
    END IF
    IF (ALLOCATED(helper%inner_recv_request_list)) THEN
       DEALLOCATE(helper%inner_recv_request_list)
    END IF
    IF (ALLOCATED(helper%data_send_request_list)) THEN
       DEALLOCATE(helper%data_send_request_list)
    END IF
    IF (ALLOCATED(helper%data_recv_request_list)) THEN
       DEALLOCATE(helper%data_recv_request_list)
    END IF
#endif

  END SUBROUTINE ReduceAndComposeMatrixCleanup_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third function to call, finishes setting up the matrices.
  PURE SUBROUTINE ReduceAndComposeMatrixCleanup_lsc(matrix, gathered_matrix, &
       & helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper

#include "comm_includes/ReduceAndComposeMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%outer_send_request_list)) THEN
       DEALLOCATE(helper%outer_send_request_list)
    END IF
    IF (ALLOCATED(helper%outer_recv_request_list)) THEN
       DEALLOCATE(helper%outer_recv_request_list)
    END IF
    IF (ALLOCATED(helper%inner_send_request_list)) THEN
       DEALLOCATE(helper%inner_send_request_list)
    END IF
    IF (ALLOCATED(helper%inner_recv_request_list)) THEN
       DEALLOCATE(helper%inner_recv_request_list)
    END IF
    IF (ALLOCATED(helper%data_send_request_list)) THEN
       DEALLOCATE(helper%data_send_request_list)
    END IF
    IF (ALLOCATED(helper%data_recv_request_list)) THEN
       DEALLOCATE(helper%data_recv_request_list)
    END IF
#endif

  END SUBROUTINE ReduceAndComposeMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Reduce and sum the matrices in one step. If you use this method, you
  !> lose the opportunity for overlapping communication.
  SUBROUTINE ReduceAndComposeMatrix_lsr(matrix, comm, gathered_matrix)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)          :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT) :: gathered_matrix
    !! Local Variables
    TYPE(ReduceHelper_t) :: helper

#include "comm_includes/ReduceAndComposeMatrix.f90"

  END SUBROUTINE ReduceAndComposeMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Reduce and sum the matrices in one step. If you use this method, you
  !> lose the opportunity for overlapping communication.
  SUBROUTINE ReduceAndComposeMatrix_lsc(matrix, comm, gathered_matrix)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)          :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix
    !! Local Variables
    TYPE(ReduceHelper_t) :: helper

#include "comm_includes/ReduceAndComposeMatrix.f90"

  END SUBROUTINE ReduceAndComposeMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceAndSumMatrixSizes_lsr(matrix, comm, gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The  helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndSumMatrixSizes_sendrecv.f90"
#else
#include "comm_includes/ReduceAndSumMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceAndSumMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceAndSumMatrixSizes_lsc(matrix, comm, gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndSumMatrixSizes_sendrecv.f90"
#else
#include "comm_includes/ReduceAndSumMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceAndSumMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  SUBROUTINE ReduceAndSumMatrixData_lsr(matrix, comm, gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndSumMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTREAL, &
            & II - 1, 4, comm, helper%data_send_request_list(II), ierr)
       istart = helper%displacement(II) + 1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, MPINTREAL, &
            & II - 1, 4, comm, helper%data_recv_request_list(II), ierr)
    END DO
#else
#include "comm_includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, comm, helper%data_request, ierr)
#endif
  END SUBROUTINE ReduceAndSumMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  SUBROUTINE ReduceAndSumMatrixData_lsc(matrix, comm, gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
#include "comm_includes/ReduceAndSumMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTCOMPLEX, &
            & II - 1, 4, comm, helper%data_send_request_list(II), ierr)
       istart = helper%displacement(II) + 1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, &
            & MPINTCOMPLEX, II - 1, 4, comm, &
            & helper%data_recv_request_list(II), ierr)
    END DO
#else
#include "comm_includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, comm, helper%data_request, ierr)
#endif
  END SUBROUTINE ReduceAndSumMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finally routine to sum up the matrices.
  PURE SUBROUTINE ReduceAndSumMatrixCleanup_lsr(matrix, gathered_matrix, &
       & threshold, helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The gathered_matrix the matrix being gathered.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The threshold the threshold for flushing values.
    REAL(NTREAL), INTENT(IN)            :: threshold
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_lsr) :: acc_matrix, sum_matrix

#include "comm_includes/ReduceAndSumMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%outer_send_request_list)) THEN
       DEALLOCATE(helper%outer_send_request_list)
    END IF
    IF (ALLOCATED(helper%outer_recv_request_list)) THEN
       DEALLOCATE(helper%outer_recv_request_list)
    END IF
    IF (ALLOCATED(helper%inner_send_request_list)) THEN
       DEALLOCATE(helper%inner_send_request_list)
    END IF
    IF (ALLOCATED(helper%inner_recv_request_list)) THEN
       DEALLOCATE(helper%inner_recv_request_list)
    END IF
    IF (ALLOCATED(helper%data_send_request_list)) THEN
       DEALLOCATE(helper%data_send_request_list)
    END IF
    IF (ALLOCATED(helper%data_recv_request_list)) THEN
       DEALLOCATE(helper%data_recv_request_list)
    END IF
#endif
  END SUBROUTINE ReduceAndSumMatrixCleanup_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finally routine to sum up the matrices.
  PURE SUBROUTINE ReduceAndSumMatrixCleanup_lsc(matrix, gathered_matrix, &
       & threshold, helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !> The threshold the threshold for flushing values.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !> The threshold the threshold for flushing values.
    REAL(NTREAL), INTENT(IN)            :: threshold
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_lsc) :: acc_matrix, sum_matrix

#include "comm_includes/ReduceAndSumMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%outer_send_request_list)) THEN
       DEALLOCATE(helper%outer_send_request_list)
    END IF
    IF (ALLOCATED(helper%outer_recv_request_list)) THEN
       DEALLOCATE(helper%outer_recv_request_list)
    END IF
    IF (ALLOCATED(helper%inner_send_request_list)) THEN
       DEALLOCATE(helper%inner_send_request_list)
    END IF
    IF (ALLOCATED(helper%inner_recv_request_list)) THEN
       DEALLOCATE(helper%inner_recv_request_list)
    END IF
    IF (ALLOCATED(helper%data_send_request_list)) THEN
       DEALLOCATE(helper%data_send_request_list)
    END IF
    IF (ALLOCATED(helper%data_recv_request_list)) THEN
       DEALLOCATE(helper%data_recv_request_list)
    END IF
#endif
  END SUBROUTINE ReduceAndSumMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Reduce and sum the matrices in one step. If you use this method, you
  !> lose the opportunity for overlapping communication.
  SUBROUTINE ReduceAndSumMatrix_lsr(matrix, comm, gathered_matrix, threshold)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The gathered_matrix the matrix being gathered.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The threshold the threshold for flushing values.
    REAL(NTREAL), INTENT(IN)            :: threshold
    !! Local Data
    TYPE(ReduceHelper_t) :: helper

#include "comm_includes/ReduceAndSumMatrix.f90"
  END SUBROUTINE ReduceAndSumMatrix_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Reduce and sum the matrices in one step. If you use this method, you
  !> lose the opportunity for overlapping communication.
  SUBROUTINE ReduceAndSumMatrix_lsc(matrix, comm, gathered_matrix, threshold)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: comm
    !> The threshold the threshold for flushing values.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !> The threshold the threshold for flushing values.
    REAL(NTREAL), INTENT(IN)            :: threshold
    !! Local Data
    TYPE(ReduceHelper_t) :: helper

#include "comm_includes/ReduceAndSumMatrix.f90"
  END SUBROUTINE ReduceAndSumMatrix_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  FUNCTION TestReduceSizeRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    LOGICAL :: send_request_completed, recv_request_completed
    CALL MPI_Testall(SIZE(helper%outer_send_request_list), &
         & helper%outer_send_request_list, send_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    CALL MPI_Testall(SIZE(helper%outer_recv_request_list), &
         & helper%outer_recv_request_list, recv_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    request_completed = send_request_completed .AND. recv_request_completed
#else
    CALL MPI_Test(helper%outer_request, request_completed, &
         & MPI_STATUS_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  FUNCTION TestReduceInnerRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    LOGICAL :: send_request_completed, recv_request_completed
    CALL MPI_Testall(SIZE(helper%inner_send_request_list), &
         & helper%inner_send_request_list, send_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    CALL MPI_Testall(SIZE(helper%inner_recv_request_list), &
         & helper%inner_recv_request_list, recv_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    request_completed = send_request_completed .AND. recv_request_completed
#else
    CALL MPI_Test(helper%inner_request, request_completed, &
         & MPI_STATUS_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceInnerRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  FUNCTION TestReduceDataRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    LOGICAL :: send_request_completed, recv_request_completed
    CALL MPI_Testall(SIZE(helper%data_send_request_list), &
         & helper%data_send_request_list, send_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    CALL MPI_Testall(SIZE(helper%data_recv_request_list), &
         & helper%data_recv_request_list, recv_request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
    request_completed = send_request_completed .AND. recv_request_completed
#else
    CALL MPI_Test(helper%data_request, request_completed, &
         & MPI_STATUS_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceDataRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixReduceModule
