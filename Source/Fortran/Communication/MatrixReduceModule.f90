!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for reducing matrices across processes.
MODULE MatrixReduceModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE MatrixSRAlgebraModule, ONLY : IncrementMatrixS
  USE MatrixSRModule, ONLY : Matrix_sr, ConstructEmptyMatrixS, &
       & DestructMatrixS
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
     !! These are for summing only
     !> A buffer for gathering values.
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: value_buffer
     !> A buffer for gathering inner indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index_buffer
     !> A buffer for gathering outer indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index_buffer
  END TYPE ReduceHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndComposeData
  PUBLIC :: ReduceAndComposeCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ReduceAndSumData
  PUBLIC :: ReduceAndSumCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestReduceSizeRequest
  PUBLIC :: TestReduceOuterRequest
  PUBLIC :: TestReduceInnerRequest
  PUBLIC :: TestReduceDataRequest
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceSizes(matrix, communicator, helper)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN)      :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: grid_error

    CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

    !! Gather Information About Other Processes
    ALLOCATE(helper%values_per_process(helper%comm_size))
    CALL MPI_IAllGather(SIZE(matrix%values),1,MPI_INT,&
         & helper%values_per_process,1,MPI_INT,communicator, &
         & helper%size_request, grid_error)
  END SUBROUTINE ReduceSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second function to call, will gather the data and align it one matrix
  !! next to another.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] gathered_matrix the matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndComposeData(matrix,communicator,gathered_matrix,helper)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN)      :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(Matrix_sr), INTENT(INOUT)   :: gathered_matrix
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: grid_error
    INTEGER :: counter
    INTEGER :: total_values

    !! Build Displacement List
    ALLOCATE(helper%displacement(helper%comm_size))
    helper%displacement(1) = 0
    DO counter = 2, SIZE(helper%displacement)
       helper%displacement(counter) = helper%displacement(counter-1) + &
            & helper%values_per_process(counter-1)
    END DO

    !! Build Storage
    CALL ConstructEmptySparseMatrix(gathered_matrix, &
         & matrix%columns*helper%comm_size,matrix%rows)
    total_values = SUM(helper%values_per_process)
    ALLOCATE(gathered_matrix%values(total_values))
    ALLOCATE(gathered_matrix%inner_index(total_values))
    gathered_matrix%outer_index(1) = 0

    !! MPI Calls
    CALL MPI_IAllGatherv(matrix%values,SIZE(matrix%values),MPINTREAL,&
         & gathered_matrix%values, helper%values_per_process, &
         & helper%displacement, MPINTREAL, communicator, &
         & helper%data_request, grid_error)
    CALL MPI_IAllGatherv(matrix%inner_index,SIZE(matrix%values),MPI_INT, &
         & gathered_matrix%inner_index, helper%values_per_process, &
         & helper%displacement, MPI_INT, communicator, &
         & helper%inner_request, grid_error)
    CALL MPI_IAllGather(matrix%outer_index(2:), matrix%columns,&
         & MPI_INT, gathered_matrix%outer_index(2:), &
         & matrix%columns, MPI_INT, communicator, helper%outer_request, &
         & grid_error)
  END SUBROUTINE ReduceAndComposeData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third function to call, finishes setting up the matrices.
  !! @param[in] matrix to send.
  !! @param[in] gathered_matrix matrix we are gathering.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndComposeCleanup(matrix, gathered_matrix, helper)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN)      :: matrix
    TYPE(Matrix_sr), INTENT(INOUT)   :: gathered_matrix
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: counter, inner_counter
    INTEGER :: temp_offset

    !! Sum Up The Outer Indices
    DO counter = 1, helper%comm_size - 1
       temp_offset = counter*matrix%columns+1
       DO inner_counter = 1, matrix%columns
          gathered_matrix%outer_index(temp_offset+inner_counter) = &
               & gathered_matrix%outer_index(temp_offset) + &
               & gathered_matrix%outer_index(temp_offset+inner_counter)
       END DO
    END DO
    DEALLOCATE(helper%values_per_process)
    DEALLOCATE(helper%displacement)
  END SUBROUTINE ReduceAndComposeCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine to call for gathering and summing up the data.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE ReduceAndSumData(matrix,communicator,helper)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN)      :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: grid_error
    INTEGER :: counter
    INTEGER :: sum_total_values, sum_outer_indices

    ALLOCATE(helper%displacement(helper%comm_size))

    !! Build Displacement List
    helper%displacement(1) = 0
    DO counter = 2, SIZE(helper%displacement)
       helper%displacement(counter) = helper%displacement(counter-1) + &
            & helper%values_per_process(counter-1)
    END DO

    !! Build Storage
    sum_total_values = SUM(helper%values_per_process)
    sum_outer_indices = (matrix%columns+1)*helper%comm_size
    ALLOCATE(helper%value_buffer(sum_total_values))
    ALLOCATE(helper%inner_index_buffer(sum_total_values))
    ALLOCATE(helper%outer_index_buffer(sum_outer_indices+1))

    !! MPI Calls
    CALL MPI_IAllGatherv(matrix%values,SIZE(matrix%values),MPINTREAL,&
         & helper%value_buffer, helper%values_per_process, helper%displacement,&
         & MPINTREAL, communicator, helper%data_request, grid_error)
    CALL MPI_IAllGatherv(matrix%inner_index,SIZE(matrix%values),MPI_INT, &
         & helper%inner_index_buffer, helper%values_per_process, &
         & helper%displacement, MPI_INT, communicator, helper%inner_request, &
         & grid_error)
    CALL MPI_IAllGather(matrix%outer_index, matrix%columns+1,&
         & MPI_INT, helper%outer_index_buffer, matrix%columns+1, MPI_INT, &
         & communicator, helper%outer_request, grid_error)
  END SUBROUTINE ReduceAndSumData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finally routine to sum up the matrices.
  !! @param[in] matrix to send.
  !! @param[inout] gathered_matrix the matrix being gathered.
  !! @param[in] threshold the threshold for flushing values.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE ReduceAndSumCleanup(matrix,gathered_matrix, threshold, helper)
    !! Parameters
    TYPE(Matrix_sr), INTENT(IN)      :: matrix
    TYPE(Matrix_sr), INTENT(INOUT)   :: gathered_matrix
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(ReduceHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_sr) :: temporary_matrix
    INTEGER :: counter
    INTEGER :: temporary_total_values

    !! Build Matrix Objects
    CALL ConstructEmptyMatrixS(temporary_matrix,matrix%columns,matrix%rows)
    CALL ConstructEmptyMatrixS(gathered_matrix,matrix%columns,matrix%rows)

    !! Sum
    DO counter = 1, helper%comm_size
       temporary_total_values = helper%values_per_process(counter)
       ALLOCATE(temporary_matrix%values(temporary_total_values))
       ALLOCATE(temporary_matrix%inner_index(temporary_total_values))
       temporary_matrix%values = helper%value_buffer( &
            & helper%displacement(counter)+1: &
            & helper%displacement(counter) + helper%values_per_process(counter))
       temporary_matrix%inner_index = helper%inner_index_buffer( &
            & helper%displacement(counter)+1: &
            & helper%displacement(counter) + helper%values_per_process(counter))
       temporary_matrix%outer_index = helper%outer_index_buffer(&
            & (matrix%columns+1)*(counter-1)+1:(matrix%columns+1)*(counter))
       IF (counter .EQ. helper%comm_size) THEN
          CALL IncrementMatrixS(temporary_matrix,gathered_matrix,&
               & threshold_in=threshold)
       ELSE
          CALL IncrementMatrixS(temporary_matrix,gathered_matrix,&
               & threshold_in=REAL(0.0,NTREAL))
       END IF
       DEALLOCATE(temporary_matrix%values)
       DEALLOCATE(temporary_matrix%inner_index)
    END DO
    CALL DestructMatrixS(temporary_matrix)
    DEALLOCATE(helper%value_buffer)
    DEALLOCATE(helper%inner_index_buffer)
    DEALLOCATE(helper%outer_index_buffer)
    DEALLOCATE(helper%values_per_process)
    DEALLOCATE(helper%displacement)
  END SUBROUTINE ReduceAndSumCleanup
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
