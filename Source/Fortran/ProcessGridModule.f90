!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to manage the process grid.
MODULE ProcessGridModule
  USE DataTypesModule, ONLY : MPITypeInfoInit
  USE LoggingModule
  USE MPI
  USE ISO_C_BINDING, ONLY : c_int, c_bool
#ifdef _OPENMP
  USE omp_lib, ONLY : omp_get_num_threads
#endif
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, PUBLIC :: ProcessGrid_t
     !! Describe the grid
     INTEGER, PUBLIC :: total_processors !< total processors in the grid.
     INTEGER, PUBLIC :: num_process_rows !< number of rows in the grid.
     INTEGER, PUBLIC :: num_process_columns !< number of columns in the grid.
     INTEGER, PUBLIC :: num_process_slices !< number of 2D slices in the grid.
     INTEGER, PUBLIC :: slice_size !< the size of a 2D slice.
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
  END TYPE ProcessGrid_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE(ProcessGrid_t), PUBLIC :: global_grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructProcessGrid
  PUBLIC :: IsRoot
  PUBLIC :: SplitProcessGrid
  !! Accessors for grid information
  PUBLIC :: GetMySlice
  PUBLIC :: GetMyRow
  PUBLIC :: GetMyColumn
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup the default process grid.
  !! @param[in] world_comm_ a communicator that every process in the grid is
  !! a part of.
  !! @param[in] process_rows_ number of grid rows.
  !! @param[in] process_columns_ number of grid columns.
  !! @param[in] process_slices_ number of grid slices.
  !! @param[in] be_verbose_in set true to print process grid info.
  SUBROUTINE ConstructProcessGrid(world_comm_, process_rows_, &
       & process_columns_, process_slices_, be_verbose_in)
    !! Parameters
    INTEGER(kind=c_int), INTENT(IN) :: world_comm_
    INTEGER(kind=c_int), INTENT(IN) :: process_rows_
    INTEGER(kind=c_int), INTENT(IN) :: process_columns_
    INTEGER(kind=c_int), INTENT(IN) :: process_slices_
    LOGICAL(kind=c_bool), INTENT(IN), OPTIONAL :: be_verbose_in
    LOGICAL :: be_verbose

    !! Process Optional Parameters
    IF (PRESENT(be_verbose_in)) THEN
       be_verbose = be_verbose_in
    ELSE
       be_verbose = .FALSE.
    END IF

    CALL ConstructNewProcessGrid(global_grid, world_comm_, process_rows_, &
         & process_columns_, process_slices_)

    CALL MPITypeInfoInit()

    !! Report
    IF (IsRoot(global_grid)) THEN
       CALL ActivateLogger
    END IF
    IF (be_verbose) THEN
       CALL WriteHeader("Process Grid")
       CALL EnterSubLog
       CALL WriteListElement(key = "Process Rows", &
            & int_value_in=global_grid%num_process_rows)
       CALL WriteListElement(key = "Process Columns", &
            & int_value_in=global_grid%num_process_columns)
       CALL WriteListElement(key = "Process Slices", &
            & int_value_in=global_grid%num_process_slices)
       CALL WriteListElement(key = "Column Blocks", &
            & int_value_in=global_grid%number_of_blocks_columns)
       CALL WriteListElement(key = "Row Blocks", &
            & int_value_in=global_grid%number_of_blocks_rows)
       CALL ExitSubLog
    END IF

  END SUBROUTINE ConstructProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a process grid.
  !! @param[out] grid to construct
  !! @param[in] world_comm_ a communicator that every process in the grid is
  !! a part of.
  !! @param[in] process_rows_ number of grid rows.
  !! @param[in] process_columns_ number of grid columns.
  !! @param[in] process_slices_ number of grid slices.
  SUBROUTINE ConstructNewProcessGrid(grid, world_comm_, process_rows_, &
       & process_columns_, process_slices_)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(INOUT) :: grid
    INTEGER(kind=c_int), INTENT(IN) :: world_comm_
    INTEGER(kind=c_int), INTENT(IN) :: process_rows_
    INTEGER(kind=c_int), INTENT(IN) :: process_columns_
    INTEGER(kind=c_int), INTENT(IN) :: process_slices_
    INTEGER :: ierr
    !! Local Data
    INTEGER :: column_block_multiplier
    INTEGER :: row_block_multiplier
    INTEGER :: row_counter, column_counter
    INTEGER :: num_threads

    CALL MPI_COMM_DUP(world_comm_, grid%global_comm, ierr)
    !! Grid Dimensions
    grid%num_process_rows = process_rows_
    grid%num_process_columns = process_columns_
    grid%num_process_slices = process_slices_
    CALL MPI_COMM_SIZE(grid%global_comm, grid%total_processors, ierr)
    grid%slice_size = grid%total_processors/grid%num_process_slices

    !! Grid ID
    CALL MPI_COMM_RANK(grid%global_comm,grid%global_rank,ierr)
    grid%my_slice = grid%global_rank/grid%slice_size
    grid%my_row = MOD(grid%global_rank, grid%slice_size)/grid%num_process_columns
    grid%my_column = MOD(grid%global_rank, grid%num_process_columns)

    !! Grid Communicators
    CALL MPI_COMM_SPLIT(grid%global_comm, grid%my_slice, grid%global_rank, &
         grid%within_slice_comm, ierr)
    CALL MPI_COMM_RANK(grid%within_slice_comm, grid%within_slice_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%global_comm, grid%within_slice_rank, &
         & grid%global_rank, grid%between_slice_comm, ierr)
    CALL MPI_COMM_RANK(grid%between_slice_comm, grid%between_slice_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_row, grid%global_rank, &
         & grid%row_comm, ierr)
    CALL MPI_COMM_RANK(grid%row_comm, grid%row_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_column, &
         & grid%global_rank, grid%column_comm, ierr)
    CALL MPI_COMM_RANK(grid%column_comm, grid%column_rank, ierr)

    !! Blocking Information
    column_block_multiplier = (grid%num_process_rows/grid%num_process_columns)*&
         & grid%num_process_slices
    IF (column_block_multiplier .EQ. 0) THEN
       column_block_multiplier = 1*grid%num_process_slices
    END IF
    row_block_multiplier = (grid%num_process_columns/grid%num_process_rows)* &
         & grid%num_process_slices
    IF (row_block_multiplier .EQ. 0) THEN
       row_block_multiplier = 1*grid%num_process_slices
    END IF

    !! The rule right now seems to be to have at least half as many blocks as
    !! threads.
#if defined NOBLOCK
    grid%block_multiplier = 1
#elif defined _OPENMP
    !$omp PARALLEL
    num_threads = omp_get_num_threads()
    !$omp end PARALLEL
    grid%block_multiplier = num_threads/&
         & (column_block_multiplier+row_block_multiplier)
    IF (grid%block_multiplier .EQ. 0) THEN
       grid%block_multiplier = 1
    END IF
#else
    grid%block_multiplier = 1
#endif

    grid%number_of_blocks_columns = &
         & column_block_multiplier*grid%block_multiplier
    grid%number_of_blocks_rows = &
         & row_block_multiplier*grid%block_multiplier

    !! Create Blocked Communicators
    ALLOCATE(grid%blocked_row_comm(grid%number_of_blocks_rows))
    ALLOCATE(grid%blocked_column_comm(grid%number_of_blocks_columns))
    ALLOCATE(grid%blocked_within_slice_comm(grid%number_of_blocks_columns, &
         & grid%number_of_blocks_rows))
    ALLOCATE(grid%blocked_between_slice_comm(grid%number_of_blocks_columns, &
         & grid%number_of_blocks_rows))

    DO row_counter=1,grid%number_of_blocks_rows
       DO column_counter=1,grid%number_of_blocks_columns
          CALL MPI_COMM_SPLIT(grid%global_comm, grid%my_slice, &
               & grid%global_rank, &
               & grid%blocked_within_slice_comm(column_counter,row_counter), &
               & ierr)
          CALL MPI_COMM_SPLIT(grid%global_comm, grid%within_slice_rank, &
               & grid%global_rank, &
               & grid%blocked_between_slice_comm(column_counter,row_counter), &
               & ierr)
       END DO
    END DO
    DO column_counter=1,grid%number_of_blocks_columns
       CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_column, &
            & grid%global_rank, &
            & grid%blocked_column_comm(column_counter), ierr)
    END DO
    DO row_counter=1,grid%number_of_blocks_rows
       CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_row, &
            & grid%global_rank, grid%blocked_row_comm(row_counter), ierr)
    END DO

  END SUBROUTINE ConstructNewProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a process grid, this splits it into two grids of even size
  SUBROUTINE SplitProcessGrid(old_grid, new_grid, my_color, split_slice, &
    & between_grid_comm)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(INOUT) :: old_grid
    TYPE(ProcessGrid_t), INTENT(INOUT) :: new_grid
    INTEGER, INTENT(OUT) :: my_color
    LOGICAL, INTENT(OUT) :: split_slice
    INTEGER, INTENT(OUT) :: between_grid_comm
    !! Local Variables - new grid
    INTEGER :: new_comm
    INTEGER :: rows, cols, slices
    INTEGER :: midpoint
    !! For Between Comm
    INTEGER :: between_color, between_rank
    INTEGER :: left_grid_size
    !! Temporary
    INTEGER :: ierr

    split_slice = .FALSE.
    !! Handle base case
    IF (old_grid%total_processors .EQ. 1) THEN
       rows = 1
       cols = 1
       slices = 1
       my_color = 0
       between_rank = 0
    !! First preferentially try to split along slices
    ELSE IF (old_grid%num_process_slices .GT. 1) THEN
       midpoint = old_grid%num_process_slices/2
       cols = old_grid%num_process_columns
       rows = old_grid%num_process_rows
       IF (old_grid%my_slice .LT. midpoint) THEN
          my_color = 0
          slices = midpoint
       ELSE
          my_color = 1
          slices = old_grid%num_process_slices - midpoint
       END IF
       between_rank = old_grid%my_slice
       split_slice = .TRUE.
       left_grid_size = midpoint*cols*rows
    !! Next try to split the bigger direction
    ELSE IF (old_grid%num_process_rows .GT. old_grid%num_process_columns) THEN
       midpoint = old_grid%num_process_rows/2
       cols = old_grid%num_process_columns
       slices = 1
       IF (old_grid%my_row .LT. midpoint) THEN
          my_color = 0
          rows = midpoint
       ELSE
          my_color = 1
          rows = old_grid%num_process_rows - midpoint
       END IF
       between_rank = old_grid%my_row
       left_grid_size = midpoint*cols*slices
    !! Default Case
    ELSE
       midpoint = old_grid%num_process_columns/2
       slices = 1
       rows = old_grid%num_process_rows
       IF (old_grid%my_column .LT. midpoint) THEN
          my_color = 0
          cols = midpoint
       ELSE
          my_color = 1
          cols = old_grid%num_process_columns - midpoint
       END IF
       between_rank = old_grid%my_column
       left_grid_size = midpoint*slices*rows
    END IF

    !! Construct
    CALL MPI_COMM_SPLIT(old_grid%global_comm, my_color, old_grid%global_rank, &
         & new_comm, ierr)
    CALL ConstructNewProcessGrid(new_grid, new_comm, rows, cols, slices)

    !! For sending data between grids
    between_color = MOD(new_grid%global_rank, left_grid_size)
    CALL MPI_COMM_SPLIT(old_grid%global_comm, between_color, between_rank, &
         & between_grid_comm, ierr)

  END SUBROUTINE SplitProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Check if the current process is the root process.
  !! @return true if the current process is root.
  FUNCTION IsRoot(grid) RESULT(is_root)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(IN) :: grid
    LOGICAL :: is_root
    IF (grid%global_rank == 0) THEN
       is_root = .TRUE.
    ELSE
       is_root = .FALSE.
    END IF
  END FUNCTION IsRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the slice of the current process.
  !! @param[in] grid the process grid (optional).
  !! @return slice number of the current process.
  FUNCTION GetMySlice(grid) RESULT(return_val)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid

    INTEGER :: return_val
    IF (PRESENT(grid)) THEN
       return_val = grid%my_slice
    ELSE
       return_val = global_grid%my_slice
    END IF
  END FUNCTION GetMySlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the column of the current process.
  !! @param[in] grid the process grid (optional).
  !! @return column number of the current process.
  FUNCTION GetMyColumn(grid) RESULT(return_val)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid
    INTEGER :: return_val

    IF (PRESENT(grid)) THEN
       return_val = grid%my_column
    ELSE
       return_val = global_grid%my_column
    END IF
  END FUNCTION GetMyColumn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the row of the current process.
  !! @param[in] grid the process grid (optional).
  !! @return row number of the current process.
  FUNCTION GetMyRow(grid) RESULT(return_val)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid
    INTEGER :: return_val

    IF (PRESENT(grid)) THEN
       return_val = grid%my_row
    ELSE
       return_val = global_grid%my_row
    END IF
  END FUNCTION GetMyRow
END MODULE ProcessGridModule
