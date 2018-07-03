!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for reducing matrices across processes.
MODULE MatrixReduceModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE MatrixSAlgebraModule, ONLY : IncrementMatrix
  USE MatrixSModule, ONLY : Matrix_lsr, Matrix_lsc, ConstructEmptyMatrix, &
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
     !> A status object for gathering the sizes.
     INTEGER :: size_status(MPI_STATUS_SIZE)
     !> A request object for gathering outer indices.
     INTEGER :: outer_request
     !> A status object for gathering outer indices.
     INTEGER :: outer_status(MPI_STATUS_SIZE)
     !> A request object for gathering inner indices.
     INTEGER :: inner_request
     !> A status object for gathering inner indices.
     INTEGER :: inner_status(MPI_STATUS_SIZE)
     !> A request object for gathering data.
     INTEGER :: data_request
     !> A status object for gathering data.
     INTEGER :: data_status(MPI_STATUS_SIZE)
     !> The error code after an MPI call.
     INTEGER :: error_code
     !> Number of values to gather from each process.
     INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_process
     !> The displacements for where those gathered values should go.
     INTEGER, DIMENSION(:), ALLOCATABLE :: displacement
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
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceMatrixSizes_lsr(matrix, communicator, helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)      :: matrix

    INCLUDE "includes/ReduceMatrixSizes.f90"
  END SUBROUTINE ReduceMatrixSizes_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceMatrixSizes_lsc(matrix, communicator, helper)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)      :: matrix

    INCLUDE "includes/ReduceMatrixSizes.f90"
  END SUBROUTINE ReduceMatrixSizes_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align it one matrix
  !! next to another.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] gathered_matrix the matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndComposeMatrixData_lsr(matrix, communicator, &
       & gathered_matrix, helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    TYPE(Matrix_lsr), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, communicator, helper%data_request, &
         & grid_error)
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
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndComposeMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, communicator, helper%data_request, &
         & grid_error)
  END SUBROUTINE ReduceAndComposeMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third function to call, finishes setting up the matrices.
  !! @param[in] matrix to send.
  !! @param[in] gathered_matrix matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndComposeMatrixCleanup_lsr(matrix, gathered_matrix, &
       & helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    TYPE(Matrix_lsr), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndComposeMatrixCleanup.f90"

  END SUBROUTINE ReduceAndComposeMatrixCleanup_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third function to call, finishes setting up the matrices.
  !! @param[in] matrix to send.
  !! @param[in] gathered_matrix matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndComposeMatrixCleanup_lsc(matrix, gathered_matrix, &
       & helper)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndComposeMatrixCleanup.f90"

  END SUBROUTINE ReduceAndComposeMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndSumMatrixData_lsr(matrix, gathered_matrix, communicator, &
       & helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    TYPE(Matrix_lsr), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, communicator, helper%data_request, &
         & grid_error)
  END SUBROUTINE ReduceAndSumMatrixData_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndSumMatrixData_lsc(matrix, gathered_matrix, communicator, &
       & helper)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)    :: matrix
    TYPE(Matrix_lsc), INTENT(INOUT) :: gathered_matrix

    INCLUDE "includes/ReduceAndSumMatrixData.f90"
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTCOMPLEX,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTCOMPLEX, communicator, &
         & helper%data_request, grid_error)
  END SUBROUTINE ReduceAndSumMatrixData_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finally routine to sum up the matrices.
  !! @param[in] matrix to send.
  !! @param[inout] gathered_matrix the matrix being gathered.
  !! @param[in] threshold the threshold for flushing values.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndSumMatrixCleanup_lsr(matrix, gathered_matrix, &
       & threshold, helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)      :: matrix
    TYPE(Matrix_lsr), INTENT(INOUT)   :: gathered_matrix
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_lsr) :: temporary_matrix, sum_matrix

    INCLUDE "includes/ReduceAndSumMatrixCleanup.f90"
  END SUBROUTINE ReduceAndSumMatrixCleanup_lsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finally routine to sum up the matrices.
  !! @param[in] matrix to send.
  !! @param[inout] gathered_matrix the matrix being gathered.
  !! @param[in] threshold the threshold for flushing values.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndSumMatrixCleanup_lsc(matrix, gathered_matrix, &
       & threshold, helper)
    !! Parameters
    TYPE(Matrix_lsc), INTENT(IN)      :: matrix
    TYPE(Matrix_lsc), INTENT(INOUT)   :: gathered_matrix
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_lsc) :: temporary_matrix, sum_matrix

    INCLUDE "includes/ReduceAndSumMatrixCleanup.f90"
  END SUBROUTINE ReduceAndSumMatrixCleanup_lsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestReduceSizeRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%size_request, request_completed, &
         & helper%size_status, helper%error_code)
  END FUNCTION TestReduceSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestReduceOuterRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%outer_request, request_completed, &
         & helper%outer_status, helper%error_code)
  END FUNCTION TestReduceOuterRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestReduceInnerRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%inner_request, request_completed, &
         & helper%inner_status, helper%error_code)
  END FUNCTION TestReduceInnerRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestReduceDataRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%data_request, request_completed, &
         & helper%data_status, helper%error_code)
  END FUNCTION TestReduceDataRequest
END MODULE MatrixReduceModule
