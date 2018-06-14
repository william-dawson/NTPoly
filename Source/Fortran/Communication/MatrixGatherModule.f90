!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for gathering matrices across processes.
MODULE MatrixGatherModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE MatrixDModule, ONLY : Matrix_ldr, ConstructEmptyMatrix
  USE MatrixSModule, ONLY : Matrix_lsr, ConstructEmptyMatrix, &
       & DestructMatrix, CopyMatrix
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data structure to stores internal information about a gather call.
  TYPE, PUBLIC :: GatherHelper_t
     !> The size of this process.
     INTEGER :: comm_size
     !> The error code after an MPI call.
     INTEGER :: error_code
     !> A request object for broadcasting the sizes.
     INTEGER :: size_request
     !> A status object for gathering the sizes.
     INTEGER :: size_status(MPI_STATUS_SIZE)
     !> A request object for broadcasting outer indices.
     INTEGER :: outer_request
     !> A status object for gathering outer indices.
     INTEGER :: outer_status(MPI_STATUS_SIZE)
     !> A request object for broadcasting inner indices.
     INTEGER :: inner_request
     !> A status object for gathering inner indices.
     INTEGER :: inner_status(MPI_STATUS_SIZE)
     !> A request object for broadcasting data.
     INTEGER :: data_request
     !> A status object for gathering data.
     INTEGER :: data_status(MPI_STATUS_SIZE)
     !> Number of values to receive, rows, columns
     INTEGER, DIMENSION(:), ALLOCATABLE :: matrix_data
     !> A buffer for gathering values.
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: value_buffer
     !> A buffer for gathering inner indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index_buffer
     !> A buffer for gathering outer indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index_buffer
     !> How many outer indices for each matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_sizes
     !> How many values for each matrix.
     INTEGER, DIMENSION(:), ALLOCATABLE :: value_sizes
     !> Displacement of the outer indices.
     INTEGER, DIMENSION(:), ALLOCATABLE :: displacement_outer
     !> Displacement of the values.
     INTEGER, DIMENSION(:), ALLOCATABLE :: displacement_values
  END TYPE GatherHelper_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: BlockingSparseMatrixGather
  PUBLIC :: GatherSparseMatrixSizes
  PUBLIC :: GatherSparseMatrixAndListData
  PUBLIC :: GatherSparseMatrixAndListCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: BlockingDenseMatrixGather
  PUBLIC :: GatherDenseMatrixSizes
  PUBLIC :: GatherDenseMatrixAndListData
  PUBLIC :: GatherDenseMatrixAndListCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TestGatherSizeRequest
  PUBLIC :: TestGatherOuterRequest
  PUBLIC :: TestGatherInnerRequest
  PUBLIC :: TestGatherDataRequest
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A blocking call to the full matrix gather routines.
  !! Each processor contributes a matrix, which are gathered on all processors
  !! in a list.
  !! @param[in] matrix the matrix to gather.
  !! @param[out] matrix_list the list of matrices to gather
  !! @param[inout] communicator the communicator to gather on.
  SUBROUTINE BlockingSparseMatrixGather(matrix, matrix_list, communicator)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: matrix
    TYPE(Matrix_lsr), DIMENSION(:), INTENT(INOUT) :: matrix_list
    INTEGER, INTENT(INOUT)           :: communicator
    !! Local Data
    TYPE(GatherHelper_t) :: helper

    !! Perfomr all the steps
    CALL GatherSparseMatrixSizes(matrix, communicator, helper)
    DO WHILE (.NOT. TestGatherSizeRequest(helper))
    END DO
    CALL GatherSparseMatrixAndListData(matrix, communicator, helper)
    DO WHILE (.NOT. TestGatherOuterRequest(helper))
    END DO
    DO WHILE (.NOT. TestGatherInnerRequest(helper))
    END DO
    DO WHILE (.NOT. TestGatherDataRequest(helper))
    END DO
    CALL GatherSparseMatrixAndListCleanup(matrix, matrix_list, helper)

  END SUBROUTINE BlockingSparseMatrixGather
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE GatherSparseMatrixSizes(matrix, communicator, helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)    :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: grid_error
    INTEGER, DIMENSION(3) :: send_data

    CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

    !! Setup the matrix data array
    ALLOCATE(helper%matrix_data(3*helper%comm_size))
    send_data(1) = SIZE(matrix%values)
    send_data(2) = matrix%rows
    send_data(3) = matrix%columns

    !! Gather Information About Other Processes
    CALL MPI_IAllGather(send_data, 3, MPI_INT,&
         & helper%matrix_data, 3, MPI_INT, communicator, &
         & helper%size_request, grid_error)

  END SUBROUTINE GatherSparseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine we gather the data into a list.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE GatherSparseMatrixAndListData(matrix,communicator,helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN)      :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Temporary Values
    INTEGER :: grid_error
    INTEGER :: II

    !! Lists of Sizes
    ALLOCATE(helper%displacement_outer(helper%comm_size))
    ALLOCATE(helper%displacement_values(helper%comm_size))
    ALLOCATE(helper%value_sizes(helper%comm_size))
    ALLOCATE(helper%outer_sizes(helper%comm_size))
    DO II = 1, helper%comm_size
       helper%value_sizes(II) = helper%matrix_data((II-1)*3 + 1)
       helper%outer_sizes(II) = helper%matrix_data((II-1)*3 + 3) + 1
    END DO

    !! Build Displacement List
    helper%displacement_values(1) = 0
    helper%displacement_outer(1) = 0
    DO II = 2, helper%comm_size
       helper%displacement_values(II) = helper%displacement_values(II-1) + &
            & helper%value_sizes(II-1)
       helper%displacement_outer(II) = helper%displacement_outer(II-1) + &
            & helper%outer_sizes(II-1)
    END DO

    !! Build Storage
    ALLOCATE(helper%value_buffer(SUM(helper%value_sizes)))
    ALLOCATE(helper%inner_index_buffer(SUM(helper%value_sizes)))
    ALLOCATE(helper%outer_index_buffer(SUM(helper%outer_sizes)))

    !! MPI Calls
    CALL MPI_IAllGatherv(matrix%values, SIZE(matrix%values), MPINTREAL,&
         & helper%value_buffer, helper%value_sizes, helper%displacement_values,&
         & MPINTREAL, communicator, helper%data_request, grid_error)
    CALL MPI_IAllGatherv(matrix%inner_index, SIZE(matrix%values), MPI_INT, &
         & helper%inner_index_buffer, helper%value_sizes, &
         & helper%displacement_values, MPI_INT, communicator, &
         & helper%inner_request, grid_error)
    CALL MPI_IAllGatherv(matrix%outer_index, SIZE(matrix%outer_index), MPI_INT,&
         & helper%outer_index_buffer, helper%outer_sizes, &
         & helper%displacement_outer, MPI_INT, communicator, &
         & helper%outer_request, grid_error)

  END SUBROUTINE GatherSparseMatrixAndListData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third a routine to cleanup the list builder.
  !! @param[in] matrix to send.
  !! @param[inout] matrix_list the list of matrices being gathered.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE GatherSparseMatrixAndListCleanup(matrix, matrix_list, helper)
    !! Parameters
    TYPE(Matrix_lsr), INTENT(IN) :: matrix
    TYPE(Matrix_lsr), DIMENSION(:), INTENT(INOUT) :: matrix_list
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    TYPE(Matrix_lsr) :: temporary_matrix
    INTEGER :: II

    !! Unpack into the list
    DO II = 1, helper%comm_size
       !! Allocate Memory
       CALL ConstructEmptyMatrix(temporary_matrix, &
            & helper%matrix_data(3*(II-1)+2), helper%matrix_data(3*(II-1)+2))
       ALLOCATE(temporary_matrix%values(helper%value_sizes(II)))
       ALLOCATE(temporary_matrix%inner_index(helper%value_sizes(II)))

       !! Copy Values
       temporary_matrix%values = helper%value_buffer( &
            & helper%displacement_values(II)+1: &
            & helper%displacement_values(II) + helper%value_sizes(II))
       temporary_matrix%inner_index = helper%inner_index_buffer( &
            & helper%displacement_values(II)+1: &
            & helper%displacement_values(II) + helper%value_sizes(II))
       temporary_matrix%outer_index = helper%outer_index_buffer( &
            & helper%displacement_outer(II)+1: &
            & helper%displacement_outer(II) + helper%outer_sizes(II))

       !! Copy from temporary to the list
       CALL CopyMatrix(temporary_matrix, matrix_list(II))
       CALL DestructMatrix(temporary_matrix)
    END DO

    !! Cleanup
    DEALLOCATE(helper%value_buffer)
    DEALLOCATE(helper%inner_index_buffer)
    DEALLOCATE(helper%outer_index_buffer)
    DEALLOCATE(helper%value_sizes)
    DEALLOCATE(helper%outer_sizes)
    DEALLOCATE(helper%displacement_values)
    DEALLOCATE(helper%displacement_outer)

  END SUBROUTINE GatherSparseMatrixAndListCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A blocking call to the full matrix gather routines.
  !! Each processor contributes a matrix, which are gathered on all processors
  !! in a list.
  !! @param[in] matrix the matrix to gather.
  !! @param[out] matrix_list the list of matrices to gather
  !! @param[inout] communicator the communicator to gather on.
  SUBROUTINE BlockingDenseMatrixGather(matrix, matrix_list, communicator)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: matrix
    TYPE(Matrix_ldr), DIMENSION(:), INTENT(INOUT) :: matrix_list
    INTEGER, INTENT(INOUT)           :: communicator
    !! Local Data
    TYPE(GatherHelper_t) :: helper

    !! Perfomr all the steps
    CALL GatherDenseMatrixSizes(matrix, communicator, helper)
    DO WHILE (.NOT. TestGatherSizeRequest(helper))
    END DO
    CALL GatherDenseMatrixAndListData(matrix, communicator, helper)
    DO WHILE (.NOT. TestGatherDataRequest(helper))
    END DO
    CALL GatherDenseMatrixAndListCleanup(matrix, matrix_list, helper)

  END SUBROUTINE BlockingDenseMatrixGather
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The first routine to call, gathers the sizes of the data to be sent.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE GatherDenseMatrixSizes(matrix, communicator, helper)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN)    :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: grid_error
    INTEGER, DIMENSION(3) :: send_data

    CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

    !! Setup the matrix data array
    IF (ALLOCATED(helper%matrix_data)) DEALLOCATE(helper%matrix_data)
    ALLOCATE(helper%matrix_data(3*helper%comm_size))
    send_data(1) = matrix%rows*matrix%columns
    send_data(2) = matrix%rows
    send_data(3) = matrix%columns

    !! Gather Information About Other Processes
    CALL MPI_IAllGather(send_data, 3, MPI_INT,&
         & helper%matrix_data, 3, MPI_INT, communicator, &
         & helper%size_request, grid_error)

  END SUBROUTINE GatherDenseMatrixSizes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Second routine we gather the data into a list.
  !! @param[in] matrix to send.
  !! @param[inout] communicator to send along.
  !! @param[inout] helper a helper associated with this gather.
  SUBROUTINE GatherDenseMatrixAndListData(matrix,communicator,helper)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN)      :: matrix
    INTEGER, INTENT(INOUT)              :: communicator
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Temporary Values
    INTEGER :: grid_error
    INTEGER :: II

    !! Lists of Sizes
    ALLOCATE(helper%displacement_values(helper%comm_size))
    ALLOCATE(helper%value_sizes(helper%comm_size))
    DO II = 1, helper%comm_size
       helper%value_sizes(II) = helper%matrix_data((II-1)*3 + 1)
    END DO

    !! Build Displacement List
    helper%displacement_values(1) = 0
    DO II = 2, helper%comm_size
       helper%displacement_values(II) = helper%displacement_values(II-1) + &
            & helper%value_sizes(II-1)
    END DO

    !! Build Storage
    ALLOCATE(helper%value_buffer(SUM(helper%value_sizes)))

    !! MPI Calls
    CALL MPI_IAllGatherv(matrix%data, matrix%rows*matrix%columns, MPINTREAL,&
         & helper%value_buffer, helper%value_sizes, helper%displacement_values,&
         & MPINTREAL, communicator, helper%data_request, grid_error)

  END SUBROUTINE GatherDenseMatrixAndListData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Third a routine to cleanup the list builder.
  !! @param[in] matrix to send.
  !! @param[inout] matrix_list the list of matrices gathered.
  !! @param[inout] helper a helper associated with this gather.
  PURE SUBROUTINE GatherDenseMatrixAndListCleanup(matrix, matrix_list, helper)
    !! Parameters
    TYPE(Matrix_ldr), INTENT(IN) :: matrix
    TYPE(Matrix_ldr), DIMENSION(:), INTENT(INOUT) :: matrix_list
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    !! Local Data
    INTEGER :: II

    !! Unpack into the list
    DO II = 1, helper%comm_size
       !! Allocate Memory
       CALL ConstructEmptyMatrix(matrix_list(II), &
            & helper%matrix_data(3*(II-1)+2), helper%matrix_data(3*(II-1)+2))

       !! Copy Values
       matrix_list(II)%data = RESHAPE(helper%value_buffer( &
            & helper%displacement_values(II)+1: &
            & helper%displacement_values(II) + helper%value_sizes(II)), &
            & SHAPE=(/matrix_list(II)%rows, matrix_list(II)%columns/))
    END DO

    !! Cleanup
    DEALLOCATE(helper%value_buffer)
    DEALLOCATE(helper%value_sizes)
    DEALLOCATE(helper%displacement_values)

  END SUBROUTINE GatherDenseMatrixAndListCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the size of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestGatherSizeRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%size_request, request_completed, &
         & helper%size_status, helper%error_code)
  END FUNCTION TestGatherSizeRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the outer indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestGatherOuterRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%outer_request, request_completed, &
         & helper%outer_status, helper%error_code)
  END FUNCTION TestGatherOuterRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the inner indices of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestGatherInnerRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%inner_request, request_completed, &
         & helper%inner_status, helper%error_code)
  END FUNCTION TestGatherInnerRequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Test if a request for the data of the matrices is complete.
  !! @param[in] helper the gatherer helper structure.
  !! @return request_completed true if the request is finished.
  FUNCTION TestGatherDataRequest(helper) RESULT(request_completed)
    !! Parameters
    TYPE(GatherHelper_t), INTENT(INOUT) :: helper
    LOGICAL :: request_completed

    CALL MPI_Test(helper%data_request, request_completed, &
         & helper%data_status, helper%error_code)
  END FUNCTION TestGatherDataRequest
END MODULE MatrixGatherModule