!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for reducing matrices across processes.
MODULE MatrixReduceModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE SMatrixAlgebraModule, ONLY : IncrementMatrix
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, ConstructEmptyMatrix, &
       & DestructMatrix, CopyMatrix
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a reduce call.
  TYPE, PUBLIC :: ReduceHelper_t
     !> Number of processors involved in this gather.
     INTEGER :: comm_size
     !> A request object for gathering the sizes.
     INTEGER :: size_request
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
     !> For mpi backup plan, a list of request object for the sizes.
     INTEGER, DIMENSION(:), ALLOCATABLE :: size_request_list
     !> For mpi backup, a list of request objets for outer indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_request_list
     !> For mpi backup, a list of request objects for inner indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_request_list
     !> For mpi backup, a list of request object for data.
     INTEGER, DIMENSION(:), ALLOCATABLE :: data_request_list
#endif
  END TYPE ReduceHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndComposeMatrixData
  PUBLIC :: ReduceAndComposeMatrixCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndSumMatrixData
  PUBLIC :: ReduceAndSumMatrixCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestReduceSizeRequest
  PUBLIC :: TestReduceOuterRequest
  PUBLIC :: TestReduceInnerRequest
  PUBLIC :: TestReduceDataRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ReduceMatrixSizes
     MODULE PROCEDURE ReduceMatrixSizes_lsr
     MODULE PROCEDURE ReduceMatrixSizes_lsc
  END INTERFACE
  INTERFACE ReduceAndComposeMatrixData
     MODULE PROCEDURE ReduceAndComposeMatrixData_lsr
     MODULE PROCEDURE ReduceAndComposeMatrixData_lsc
  END INTERFACE
  INTERFACE ReduceAndComposeMatrixCleanup
     MODULE PROCEDURE ReduceAndComposeMatrixCleanup_lsr
     MODULE PROCEDURE ReduceAndComposeMatrixCleanup_lsc
  END INTERFACE
  INTERFACE ReduceAndSumMatrixData
     MODULE PROCEDURE ReduceAndSumMatrixData_lsr
     MODULE PROCEDURE ReduceAndSumMatrixData_lsc
  END INTERFACE
  INTERFACE ReduceAndSumMatrixCleanup
     MODULE PROCEDURE ReduceAndSumMatrixCleanup_lsr
     MODULE PROCEDURE ReduceAndSumMatrixCleanup_lsc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceMatrixSizes_lsr(matrix, communicator, helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
    !> The  helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper

#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceMatrixSizes_sendrecv.f90"
#else
    INCLUDE "comm_includes/ReduceMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  SUBROUTINE ReduceMatrixSizes_lsc(matrix, communicator, helper)
    !! The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !! The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
    !! The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper

#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceMatrixSizes_sendrecv.f90"
#else
    INCLUDE "comm_includes/ReduceMatrixSizes.f90"
#endif
  END SUBROUTINE ReduceMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align it one matrix
  !! @param[inout]
  SUBROUTINE ReduceAndComposeMatrixData_lsr(matrix, communicator, &
       & gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceAndComposeMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTREAL, &
            & II-1, 4, communicator, helper%data_request_list(II), grid_error)
       istart = helper%displacement(II)+1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, MPINTREAL, &
            & II-1, 4, communicator, &
            & helper%data_request_list(helper%comm_size+II), grid_error)
    END DO
#else
    INCLUDE "comm_includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, communicator, helper%data_request, &
         & grid_error)
#endif
  END SUBROUTINE ReduceAndComposeMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align it one matrix
  !! next to another.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] gathered_matrix the matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndComposeMatrixData_lsc(matrix, communicator, &
       & gathered_matrix, helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)        :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT)     :: gathered_matrix
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceAndComposeMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTCOMPLEX, &
            & II-1, 4, communicator, helper%data_request_list(II), grid_error)
       istart = helper%displacement(II)+1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, &
            & MPINTCOMPLEX, II-1, 4, communicator, &
            & helper%data_request_list(helper%comm_size+II), grid_error)
    END DO
#else
    INCLUDE "comm_includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, communicator, &
         & helper%data_request, grid_error)
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

    INCLUDE "comm_includes/ReduceAndComposeMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%size_request_list)) THEN
       DEALLOCATE(helper%size_request_list)
    END IF
    IF (ALLOCATED(helper%outer_request_list)) THEN
       DEALLOCATE(helper%outer_request_list)
    END IF
    IF (ALLOCATED(helper%inner_request_list)) THEN
       DEALLOCATE(helper%inner_request_list)
    END IF
    IF (ALLOCATED(helper%data_request_list)) THEN
       DEALLOCATE(helper%data_request_list)
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

    INCLUDE "comm_includes/ReduceAndComposeMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%size_request_list)) THEN
       DEALLOCATE(helper%size_request_list)
    END IF
    IF (ALLOCATED(helper%outer_request_list)) THEN
       DEALLOCATE(helper%outer_request_list)
    END IF
    IF (ALLOCATED(helper%inner_request_list)) THEN
       DEALLOCATE(helper%inner_request_list)
    END IF
    IF (ALLOCATED(helper%data_request_list)) THEN
       DEALLOCATE(helper%data_request_list)
    END IF
#endif

  END SUBROUTINE ReduceAndComposeMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  SUBROUTINE ReduceAndSumMatrixData_lsr(matrix, gathered_matrix, communicator, &
       & helper)
    !> The matrix to send.
    TYPE(Matrix_lsr), INTENT(IN)        :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsr), INTENT(INOUT)     :: gathered_matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceAndSumMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTREAL, &
            & II-1, 4, communicator, helper%data_request_list(II), grid_error)
       istart = helper%displacement(II)+1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, MPINTREAL, &
            & II-1, 4, communicator, &
            & helper%data_request_list(helper%comm_size+II), grid_error)
    END DO
#else
    INCLUDE "comm_includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, communicator, helper%data_request, &
         & grid_error)
#endif
  END SUBROUTINE ReduceAndSumMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndSumMatrixData_lsc(matrix, gathered_matrix, communicator, &
       & helper)
    !> The matrix to send.
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    !> The matrix we are gathering.
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix
    !> The communicator to send along.
    INTEGER, INTENT(INOUT)              :: communicator
    !> The helper associated with this gather.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
#ifdef NOIALLGATHER
    INCLUDE "comm_includes/ReduceAndSumMatrixData_sendrecv.f90"
    DO II = 1, helper%comm_size
       CALL MPI_ISend(matrix%values, SIZE(matrix%values), MPINTCOMPLEX, &
            & II-1, 4, communicator, helper%data_request_list(II), grid_error)
       istart = helper%displacement(II)+1
       isize = helper%values_per_process(II)
       iend = istart + isize - 1
       CALL MPI_Irecv(gathered_matrix%values(istart:iend), isize, &
            & MPINTCOMPLEX, II-1, 4, communicator, &
            & helper%data_request_list(helper%comm_size+II), grid_error)
    END DO
#else
    INCLUDE "comm_includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, communicator, &
         & helper%data_request, grid_error)
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
    TYPE(Matrix_lsr) :: temporary_matrix, sum_matrix

    INCLUDE "comm_includes/ReduceAndSumMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%size_request_list)) THEN
       DEALLOCATE(helper%size_request_list)
    END IF
    IF (ALLOCATED(helper%outer_request_list)) THEN
       DEALLOCATE(helper%outer_request_list)
    END IF
    IF (ALLOCATED(helper%inner_request_list)) THEN
       DEALLOCATE(helper%inner_request_list)
    END IF
    IF (ALLOCATED(helper%data_request_list)) THEN
       DEALLOCATE(helper%data_request_list)
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
    TYPE(Matrix_lsc) :: temporary_matrix, sum_matrix

    INCLUDE "comm_includes/ReduceAndSumMatrixCleanup.f90"
#ifdef NOIALLGATHER
    IF (ALLOCATED(helper%size_request_list)) THEN
       DEALLOCATE(helper%size_request_list)
    END IF
    IF (ALLOCATED(helper%outer_request_list)) THEN
       DEALLOCATE(helper%outer_request_list)
    END IF
    IF (ALLOCATED(helper%inner_request_list)) THEN
       DEALLOCATE(helper%inner_request_list)
    END IF
    IF (ALLOCATED(helper%data_request_list)) THEN
       DEALLOCATE(helper%data_request_list)
    END IF
#endif
  END SUBROUTINE ReduceAndSumMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  FUNCTION TestReduceSizeRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    CALL MPI_Testall(SIZE(helper%size_request_list), helper%size_request_list, &
         & request_completed, MPI_STATUSES_IGNORE, helper%error_code)
#else
    CALL MPI_Test(helper%size_request, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  FUNCTION TestReduceOuterRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    CALL MPI_Testall(SIZE(helper%outer_request_list), &
         & helper%outer_request_list, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
#else
    CALL MPI_Test(helper%outer_request, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceOuterRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  FUNCTION TestReduceInnerRequest(helper) RESULT(request_completed)
    !> The gatherer helper structure.
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !> True if the request is finished.
    LOGICAL :: request_completed
#ifdef NOIALLGATHER
    CALL MPI_Testall(SIZE(helper%inner_request_list), &
         & helper%inner_request_list, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
#else
    CALL MPI_Test(helper%inner_request, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
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
    CALL MPI_Testall(SIZE(helper%data_request_list), helper%data_request_list,&
         & request_completed, MPI_STATUSES_IGNORE, helper%error_code)
#else
    CALL MPI_Test(helper%data_request, request_completed, &
         & MPI_STATUSES_IGNORE, helper%error_code)
#endif
  END FUNCTION TestReduceDataRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixReduceModule
