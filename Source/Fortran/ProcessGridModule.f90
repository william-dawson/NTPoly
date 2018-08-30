!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to manage the process grid.
MODULE ProcessGridModule
  USE DataTypesModule, ONLY : MPITypeInfoInit
  USE LoggingModule, ONLY : ActivateLogger, EnterSubLog, ExitSubLog, &
       & WriteHeader, WriteListElement
  USE NTMPIModule
  USE ISO_C_BINDING, ONLY : c_int, c_bool
#ifdef _OPENMP
  USE omp_lib, ONLY : omp_get_num_threads
#endif
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Describe the grid
  INTEGER, PUBLIC :: total_processors !< total processors in the grid.
  INTEGER, PUBLIC :: num_process_rows !< number of rows in the grid.
  INTEGER, PUBLIC :: num_process_columns !< number of columns in the grid.
  INTEGER, PUBLIC :: num_process_slices !< number of 2D slices in the grid.
  INTEGER, PUBLIC :: slice_size !< the slice of a 2D slice.
  !! Identiy current process
  INTEGER, PUBLIC :: my_slice !< which slice is the current process in.
  INTEGER, PUBLIC :: my_row !< which row is the current process in.
  INTEGER, PUBLIC :: my_column !< which column is the current process in.
  !! Ranks for communication
  INTEGER, PUBLIC :: global_rank !< current process's rank amongst processes.
  INTEGER, PUBLIC :: within_slice_rank !< rank for within slice communication.
  INTEGER, PUBLIC :: between_slice_rank !< rank for between slice communication.
  INTEGER, PUBLIC :: column_rank !< rank for within column communication.
  INTEGER, PUBLIC :: row_rank !< rank for within row communication.
  !! Communicators for communication
  INTEGER, PUBLIC :: global_comm !< communicator with every other process.
  INTEGER, PUBLIC :: row_comm !< communicator within a row.
  INTEGER, PUBLIC :: column_comm !< communicator within a column.
  INTEGER, PUBLIC :: within_slice_comm !< communicator within a slice.
  INTEGER, PUBLIC :: between_slice_comm !< communicator between slices.
  INTEGER, PUBLIC :: grid_error !< stores errors from MPI calls.
  INTEGER, PUBLIC :: RootID = 0 !< Which rank is root?
  !! Blocked communication
  INTEGER, PUBLIC :: block_multiplier !< Block scaling factor.
  INTEGER, PUBLIC :: number_of_blocks_columns !< number of column blocks.
  INTEGER, PUBLIC :: number_of_blocks_rows !< number of row blocks.
  !> blocked communicator within a row.
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: blocked_row_comm
  !> blocked communicator within a column.
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: blocked_column_comm
  !> blocked communicator within a slice.
  INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: blocked_within_slice_comm
  !> blocked communicator between slices.
  INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: blocked_between_slice_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructProcessGrid
  PUBLIC :: DestructProcessGrid
  PUBLIC :: IsRoot
  !! Accessors for grid information
  PUBLIC :: GetMySlice
  PUBLIC :: GetMyRow
  PUBLIC :: GetMyColumn
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the process grid.
  !! @param[in] world_comm_ a communicator that every process in the grid is
  !! a part of.
  !! @param[in] process_rows_ number of grid rows.
  !! @param[in] process_columns_ number of grid columns.
  !! @param[in] process_slices_ number of grid slices.
  !! @param[in] be_verbose_in set true to print process grid info.
  SUBROUTINE ConstructProcessGrid(world_comm_, process_rows_, process_columns_,&
       & process_slices_, be_verbose_in)
    !! Parameters
    INTEGER(kind=c_int), INTENT(IN) :: world_comm_
    INTEGER(kind=c_int), INTENT(IN) :: process_rows_
    INTEGER(kind=c_int), INTENT(IN) :: process_columns_
    INTEGER(kind=c_int), INTENT(IN) :: process_slices_
    LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: be_verbose_in
    INTEGER :: ierr
    !! Local Data
    INTEGER :: column_block_multiplier
    INTEGER :: row_block_multiplier
    INTEGER :: row_counter, column_counter
    INTEGER :: num_threads
    LOGICAL :: be_verbose

    !! Process Optional Parameters
    IF (PRESENT(be_verbose_in)) THEN
       be_verbose = be_verbose_in
    ELSE
       be_verbose = .FALSE.
    END IF

    CALL MPI_COMM_DUP(world_comm_,global_comm, ierr)
    !! Grid Dimensions
    num_process_rows = process_rows_
    num_process_columns = process_columns_
    num_process_slices = process_slices_
    CALL MPI_COMM_SIZE(global_comm, total_processors, ierr)
    slice_size = total_processors/num_process_slices

    !! Grid ID
    CALL MPI_COMM_RANK(global_comm,global_rank,ierr)
    my_slice = global_rank/slice_size
    my_row = MOD(global_rank, slice_size)/num_process_columns
    my_column = MOD(global_rank, num_process_columns)

    !! Grid Communicators
    CALL MPI_COMM_SPLIT(global_comm, my_slice, global_rank, within_slice_comm, &
         & ierr)
    CALL MPI_COMM_RANK(within_slice_comm, within_slice_rank,ierr)
    CALL MPI_COMM_SPLIT(global_comm, within_slice_rank, global_rank, &
         & between_slice_comm, ierr)
    CALL MPI_COMM_RANK(between_slice_comm, between_slice_rank, ierr)
    CALL MPI_COMM_SPLIT(within_slice_comm, my_row, global_rank, row_comm, &
         & ierr)
    CALL MPI_COMM_RANK(row_comm, row_rank, ierr)
    CALL MPI_COMM_SPLIT(within_slice_comm, my_column, global_rank, column_comm,&
         & ierr)
    CALL MPI_COMM_RANK(column_comm, column_rank, ierr)

    !! Blocking Information
    column_block_multiplier = (num_process_rows/num_process_columns)* &
         & num_process_slices
    IF (column_block_multiplier .EQ. 0) THEN
       column_block_multiplier = 1*num_process_slices
    END IF
    row_block_multiplier = (num_process_columns/num_process_rows)* &
         & num_process_slices
    IF (row_block_multiplier .EQ. 0) THEN
       row_block_multiplier = 1*num_process_slices
    END IF

    !! The rule right now seems to be to have at least half as many blocks as
    !! threads.
#if defined NOBLOCK
    block_multiplier = 1
#elif defined _OPENMP
    !$omp PARALLEL
    num_threads = omp_get_num_threads()
    !$omp end PARALLEL
    block_multiplier = num_threads/&
         & (column_block_multiplier+row_block_multiplier)
    IF (block_multiplier .EQ. 0) THEN
       block_multiplier = 1
    END IF
#else
    block_multiplier = 1
#endif

    number_of_blocks_columns = column_block_multiplier*block_multiplier
    number_of_blocks_rows = row_block_multiplier*block_multiplier

    !! Create Blocked Communicators
    ALLOCATE(blocked_row_comm(number_of_blocks_rows))
    ALLOCATE(blocked_column_comm(number_of_blocks_columns))
    ALLOCATE(blocked_within_slice_comm(number_of_blocks_columns, &
         & number_of_blocks_rows))
    ALLOCATE(blocked_between_slice_comm(number_of_blocks_columns, &
         & number_of_blocks_rows))

    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          CALL MPI_COMM_SPLIT(global_comm, my_slice, global_rank, &
               & blocked_within_slice_comm(column_counter,row_counter), ierr)
          CALL MPI_COMM_SPLIT(global_comm, within_slice_rank, &
               & global_rank, &
               & blocked_between_slice_comm(column_counter,row_counter), ierr)
       END DO
    END DO
    DO column_counter=1,number_of_blocks_columns
       CALL MPI_COMM_SPLIT(within_slice_comm, my_column, global_rank, &
            & blocked_column_comm(column_counter), ierr)
    END DO
    DO row_counter=1,number_of_blocks_rows
       CALL MPI_COMM_SPLIT(within_slice_comm, my_row, global_rank, &
            & blocked_row_comm(row_counter), ierr)
    END DO

    CALL MPITypeInfoInit()

    !! Report
    IF (IsRoot()) THEN
       CALL ActivateLogger
    END IF
    IF (be_verbose) THEN
       CALL WriteHeader("Process Grid")
       CALL EnterSubLog
       CALL WriteListElement(key = "Process Rows", &
            & int_value_in=num_process_rows)
       CALL WriteListElement(key = "Process Columns", &
            & int_value_in=num_process_columns)
       CALL WriteListElement(key = "Process Slices", &
            & int_value_in=num_process_slices)
       CALL WriteListElement(key = "Column Blocks", &
            & int_value_in=number_of_blocks_columns)
       CALL WriteListElement(key = "Row Blocks", &
            & int_value_in=number_of_blocks_rows)
       CALL ExitSubLog
    END IF
  END SUBROUTINE ConstructProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if the current process is the root process.
  !! @return true if the current process is root.
  FUNCTION IsRoot() RESULT(is_root)
    !! Parameters
    LOGICAL :: is_root
    IF (global_rank == 0) THEN
       is_root = .TRUE.
    ELSE
       is_root = .FALSE.
    END IF
  END FUNCTION IsRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the slice of the current process.
  !! @return slice number of the current process.
  FUNCTION GetMySlice() RESULT(return_val)
    !! Parameters
    INTEGER :: return_val
    return_val = my_slice
  END FUNCTION GetMySlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the column of the current process.
  !! @return column number of the current process.
  FUNCTION GetMyColumn() RESULT(return_val)
    !! Parameters
    INTEGER :: return_val
    return_val = my_column
  END FUNCTION GetMyColumn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the row of the current process.
  !! @return row number of the current process.
  FUNCTION GetMyRow() RESULT(return_val)
    !! Parameters
    INTEGER :: return_val
    return_val = my_row
  END FUNCTION GetMyRow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct the process grid.
  SUBROUTINE DestructProcessGrid()
    !! Local Data
    INTEGER :: ierr
    INTEGER :: row_counter, column_counter

    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          CALL MPI_COMM_FREE( &
               & blocked_within_slice_comm(column_counter,row_counter), ierr)
          CALL MPI_COMM_FREE( &
               & blocked_between_slice_comm(column_counter,row_counter), ierr)
       END DO
    END DO
    DO column_counter=1,number_of_blocks_columns
       CALL MPI_COMM_FREE(blocked_column_comm(column_counter), ierr)
    END DO
    DO row_counter=1,number_of_blocks_rows
       CALL MPI_COMM_FREE(blocked_row_comm(row_counter), ierr)
    END DO

    DEALLOCATE(blocked_row_comm)
    DEALLOCATE(blocked_column_comm)
    DEALLOCATE(blocked_within_slice_comm)
    DEALLOCATE(blocked_between_slice_comm)

    CALL MPI_COMM_FREE(row_comm, ierr)
    CALL MPI_COMM_FREE(column_comm, ierr)
    CALL MPI_COMM_FREE(within_slice_comm, ierr)
    CALL MPI_COMM_FREE(between_slice_comm, ierr)
    CALL MPI_COMM_FREE(global_comm, ierr)

  END SUBROUTINE DestructProcessGrid
END MODULE ProcessGridModule
