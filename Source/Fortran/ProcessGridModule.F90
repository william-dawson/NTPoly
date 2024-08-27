!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to manage the process grid.
MODULE ProcessGridModule
  USE ErrorModule, ONLY : Error_t, ConstructError, SetGenericError
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, &
       & WriteHeader, WriteListElement
  USE NTMPIModule
#ifdef _OPENMP
  USE omp_lib, ONLY : omp_get_num_threads, omp_get_max_threads
#endif
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype which stores a process grid and all its communicators.
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
     INTEGER, PUBLIC :: global_rank !< current process rank amongst processes.
     !> rank for within slice communication.
     INTEGER, PUBLIC :: within_slice_rank
     !> rank for between slice communication.
     INTEGER, PUBLIC :: between_slice_rank
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
     !> The maximum number of openmp threads.
     INTEGER :: omp_max_threads
  END TYPE ProcessGrid_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The default process grid.
  TYPE(ProcessGrid_t), TARGET, PUBLIC, SAVE :: global_grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructProcessGrid
  PUBLIC :: ConstructNewProcessGrid
  PUBLIC :: IsRoot
  PUBLIC :: SplitProcessGrid
  PUBLIC :: CopyProcessGrid
  PUBLIC :: DestructProcessGrid
  !! Accessors for grid information
  PUBLIC :: GetMySlice
  PUBLIC :: GetMyRow
  PUBLIC :: GetMyColumn
  PUBLIC :: ComputeGridSize
  PUBLIC :: WriteProcessGridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructProcessGrid
     MODULE PROCEDURE ConstructProcessGrid_full
     MODULE PROCEDURE ConstructProcessGrid_onlyslice
  END INTERFACE ConstructProcessGrid
  INTERFACE ConstructNewProcessGrid
     MODULE PROCEDURE ConstructNewProcessGrid_full
     MODULE PROCEDURE ConstructNewProcessGrid_onlyslice
  END INTERFACE ConstructNewProcessGrid
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup the default process grid.
  SUBROUTINE ConstructProcessGrid_full(world_comm, process_rows, &
       & process_columns, process_slices)
    !> A communicator that every process in the grid is a part of.
    INTEGER, INTENT(IN) :: world_comm
    !> The number of grid rows.
    INTEGER, INTENT(IN) :: process_rows
    !> The number of grid columns.
    INTEGER, INTENT(IN) :: process_columns
    !> The number of grid slices.
    INTEGER, INTENT(IN) :: process_slices

    CALL ConstructNewProcessGrid(global_grid, world_comm, process_rows, &
         & process_columns, process_slices)
  END SUBROUTINE ConstructProcessGrid_full
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup a process grid specifying only the slices
  SUBROUTINE ConstructProcessGrid_onlyslice(world_comm, process_slices_in)
    !> A communicator that every process in the grid is a part of.
    INTEGER, INTENT(IN) :: world_comm
    !> The number of grid slices.
    INTEGER, INTENT(IN), OPTIONAL :: process_slices_in
    !! Local Data
    INTEGER :: process_rows, process_columns, process_slices
    INTEGER :: total_processors
    INTEGER :: ierr

    !! Total processors
    CALL MPI_COMM_SIZE(world_comm, total_processors, ierr)

    !! Process Optional Parameters
    IF (PRESENT(process_slices_in)) THEN
       process_slices = process_slices_in
    ELSE
       CALL ComputeNumSlices(total_processors, process_slices)
    END IF

    !! Create a 3D grid
    CALL ComputeGridSize(total_processors, process_slices, process_rows, &
         & process_columns)

    !! Now call the full setup
    CALL ConstructProcessGrid(world_comm, process_rows, process_columns, &
         & process_slices)
  END SUBROUTINE ConstructProcessGrid_onlyslice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a process grid.
  SUBROUTINE ConstructNewProcessGrid_full(grid, world_comm, process_rows, &
       & process_columns, process_slices)
    !> The grid to construct
    TYPE(ProcessGrid_t), INTENT(INOUT) :: grid
    !> A communicator that every process in the grid is a part of.
    INTEGER, INTENT(IN) :: world_comm
    !> The number of grid rows.
    INTEGER, INTENT(IN) :: process_rows
    !> The number of grid columns.
    INTEGER, INTENT(IN) :: process_columns
    !> The number of grid slices.
    INTEGER, INTENT(IN) :: process_slices
    !! Local Data
    INTEGER :: column_block_multiplier
    INTEGER :: row_block_multiplier
    INTEGER :: II, JJ
#ifdef _OPENMP
    INTEGER :: num_threads
#endif
    INTEGER :: ierr
    TYPE(Error_t) :: err

    CALL ConstructError(err)
    CALL MPI_COMM_DUP(world_comm, grid%global_comm, ierr)
    !! Grid Dimensions
    grid%num_process_rows = process_rows
    grid%num_process_columns = process_columns
    grid%num_process_slices = process_slices
    CALL MPI_COMM_SIZE(grid%global_comm, grid%total_processors, ierr)
    grid%slice_size = grid%total_processors / grid%num_process_slices

    !! Do a sanity check
    IF (grid%num_process_rows * grid%num_process_columns &
         & * grid%num_process_slices .NE. grid%total_processors) THEN
       CALL SetGenericError(err, &
            & "you did not specify a consistent process grid size", .TRUE.)
    END IF
    IF (grid%num_process_slices .GT. 1) THEN
       IF (MOD(MAX(grid%num_process_rows, grid%num_process_columns), &
            & MIN(grid%num_process_rows, grid%num_process_columns)) &
            & .NE. 0) THEN
          CALL SetGenericError(err, &
               & "if slices >1, either rows or columns must be a multiple" // &
               & "of the other.", &
               & .TRUE.)
       END IF
    END IF

    !! Grid ID
    CALL MPI_COMM_RANK(grid%global_comm,grid%global_rank, ierr)
    grid%my_slice = grid%global_rank / grid%slice_size
    grid%my_row = MOD(grid%global_rank, grid%slice_size) / &
         & grid%num_process_columns
    grid%my_column = MOD(grid%global_rank, grid%num_process_columns)

    !! Grid Communicators
    CALL MPI_COMM_SPLIT(grid%global_comm, grid%my_slice, grid%global_rank, &
         grid%within_slice_comm, ierr)
    CALL MPI_COMM_RANK(grid%within_slice_comm, grid%within_slice_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%global_comm, grid%within_slice_rank, &
         & grid%global_rank, grid%between_slice_comm, ierr)
    CALL MPI_COMM_RANK(grid%between_slice_comm, grid%between_slice_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_row, &
         & grid%global_rank, grid%row_comm, ierr)
    CALL MPI_COMM_RANK(grid%row_comm, grid%row_rank, ierr)
    CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_column, &
         & grid%global_rank, grid%column_comm, ierr)
    CALL MPI_COMM_RANK(grid%column_comm, grid%column_rank, ierr)

    !! Blocking Information
    column_block_multiplier = &
         & (grid%num_process_rows / grid%num_process_columns) * &
         & grid%num_process_slices
    IF (column_block_multiplier .EQ. 0) THEN
       column_block_multiplier = grid%num_process_slices
    END IF
    row_block_multiplier = &
         & (grid%num_process_columns / grid%num_process_rows) * &
         & grid%num_process_slices
    IF (row_block_multiplier .EQ. 0) THEN
       row_block_multiplier = 1*grid%num_process_slices
    END IF

    !! The rule right now seems to be to have at least half as many blocks as
    !! threads.
#if defined NOBLOCK
    grid%block_multiplier = 1
    grid%omp_max_threads = 1
#elif defined _OPENMP
    !$omp PARALLEL
    num_threads = omp_get_num_threads()
    grid%omp_max_threads = omp_get_max_threads()
    !$omp end PARALLEL
    grid%block_multiplier = num_threads / &
         & (column_block_multiplier + row_block_multiplier)
    IF (grid%block_multiplier .EQ. 0) THEN
       grid%block_multiplier = 1
    END IF
#else
    grid%block_multiplier = 1
#endif

    grid%number_of_blocks_columns = &
         & column_block_multiplier * grid%block_multiplier
    grid%number_of_blocks_rows = &
         & row_block_multiplier * grid%block_multiplier

    !! Create Blocked Communicators
    ALLOCATE(grid%blocked_row_comm(grid%number_of_blocks_rows))
    ALLOCATE(grid%blocked_column_comm(grid%number_of_blocks_columns))
    ALLOCATE(grid%blocked_within_slice_comm(grid%number_of_blocks_rows, &
         & grid%number_of_blocks_columns))
    ALLOCATE(grid%blocked_between_slice_comm(grid%number_of_blocks_rows, &
         & grid%number_of_blocks_columns))

    DO JJ = 1, grid%number_of_blocks_columns
       DO II = 1, grid%number_of_blocks_rows
          CALL MPI_COMM_SPLIT(grid%global_comm, grid%my_slice, &
               & grid%global_rank, grid%blocked_within_slice_comm(II, JJ), &
               & ierr)
          CALL MPI_COMM_SPLIT(grid%global_comm, grid%within_slice_rank, &
               & grid%global_rank, grid%blocked_between_slice_comm(II, JJ), &
               & ierr)
       END DO
    END DO
    DO JJ = 1, grid%number_of_blocks_columns
       CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_column, &
            & grid%global_rank, grid%blocked_column_comm(JJ), ierr)
    END DO
    DO II = 1, grid%number_of_blocks_rows
       CALL MPI_COMM_SPLIT(grid%within_slice_comm, grid%my_row, &
            & grid%global_rank, grid%blocked_row_comm(II), ierr)
    END DO

  END SUBROUTINE ConstructNewProcessGrid_full
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup a process grid specifying only the slices
  SUBROUTINE ConstructNewProcessGrid_onlyslice(grid, world_comm, &
       & process_slices_in)
    !> The grid to construct
    TYPE(ProcessGrid_t), INTENT(INOUT) :: grid
    !> A communicator that every process in the grid is a part of.
    INTEGER, INTENT(IN) :: world_comm
    !> The number of grid slices.
    INTEGER, INTENT(IN), OPTIONAL :: process_slices_in
    !! Local Data
    INTEGER :: process_rows, process_columns, process_slices
    INTEGER :: total_processors
    INTEGER :: ierr

    !! Total processors
    CALL MPI_COMM_SIZE(world_comm, total_processors, ierr)

    !! Process Optional Parameters
    IF (PRESENT(process_slices_in)) THEN
       process_slices = process_slices_in
    ELSE
       CALL ComputeNumSlices(total_processors, process_slices)
    END IF

    !! Create a 3D grid
    CALL ComputeGridSize(total_processors, process_slices, process_rows, &
         & process_columns)

    !! Now call the full setup
    CALL ConstructNewProcessGrid(grid, world_comm, process_rows, &
         & process_columns, process_slices)
  END SUBROUTINE ConstructNewProcessGrid_onlyslice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a process grid.
  !> Note that this makes a complete and independent copy of the process grid.
  !> Which of course means that whatever is currently stored in new_grid will
  !> be destroyed, so do not leave any matrices pointing to it.
  SUBROUTINE CopyProcessGrid(old_grid, new_grid)
    !> The grid to copy.
    TYPE(ProcessGrid_t), INTENT(IN) :: old_grid
    !> New_grid = old_grid
    TYPE(ProcessGrid_t), INTENT(INOUT) :: new_grid
    INTEGER :: II, JJ, ierr

    !! Safe Copy
    CALL DestructProcessGrid(new_grid)

    !! Copy all the basic data about a process
    new_grid%total_processors = old_grid%total_processors
    new_grid%num_process_rows = old_grid%num_process_rows
    new_grid%num_process_columns = old_grid%num_process_columns
    new_grid%num_process_slices = old_grid%num_process_slices
    new_grid%slice_size = old_grid%slice_size
    new_grid%my_slice = old_grid%my_slice
    new_grid%my_row = old_grid%my_row
    new_grid%my_column = old_grid%my_column
    new_grid%global_rank = old_grid%global_rank
    new_grid%within_slice_rank = old_grid%within_slice_rank
    new_grid%between_slice_rank = old_grid%between_slice_rank
    new_grid%column_rank = old_grid%column_rank
    new_grid%row_rank = old_grid%row_rank
    new_grid%block_multiplier = old_grid%block_multiplier
    new_grid%number_of_blocks_columns = old_grid%number_of_blocks_columns
    new_grid%number_of_blocks_rows = old_grid%number_of_blocks_rows
    new_grid%omp_max_threads = old_grid%omp_max_threads

    !! Allocate Blocks
    ALLOCATE(new_grid%blocked_row_comm(old_grid%number_of_blocks_rows))
    ALLOCATE(new_grid%blocked_column_comm(old_grid%number_of_blocks_columns))
    ALLOCATE(new_grid%blocked_within_slice_comm(&
         & old_grid%number_of_blocks_rows, old_grid%number_of_blocks_columns))
    ALLOCATE(new_grid%blocked_between_slice_comm( &
         & old_grid%number_of_blocks_rows, old_grid%number_of_blocks_columns))

    !! Copy the communicators
    DO II = 1, new_grid%number_of_blocks_rows
       CALL MPI_COMM_DUP(old_grid%blocked_row_comm(II), &
            & new_grid%blocked_row_comm(II), ierr)
    END DO
    DO JJ = 1, new_grid%number_of_blocks_columns
       CALL MPI_COMM_DUP(old_grid%blocked_column_comm(JJ), &
            & new_grid%blocked_column_comm(JJ), ierr)
    END DO
    DO JJ = 1, new_grid%number_of_blocks_columns
       DO II = 1,new_grid%number_of_blocks_rows
          CALL MPI_COMM_DUP(old_grid%blocked_within_slice_comm(II, JJ), &
               & new_grid%blocked_within_slice_comm(II, JJ), ierr)
       END DO
    END DO
    DO JJ = 1, new_grid%number_of_blocks_columns
       DO II = 1, new_grid%number_of_blocks_rows
          CALL MPI_COMM_DUP(old_grid%blocked_between_slice_comm(II, JJ), &
               & new_grid%blocked_between_slice_comm(II, JJ), ierr)
       END DO
    END DO

    CALL MPI_COMM_DUP(old_grid%global_comm, new_grid%global_comm, ierr)
    CALL MPI_COMM_DUP(old_grid%row_comm, new_grid%row_comm, ierr)
    CALL MPI_COMM_DUP(old_grid%column_comm, new_grid%column_comm, ierr)
    CALL MPI_COMM_DUP(old_grid%within_slice_comm, new_grid%within_slice_comm, &
         & ierr)
    CALL MPI_COMM_DUP(old_grid%between_slice_comm, &
         & new_grid%between_slice_comm, ierr)

  END SUBROUTINE CopyProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a process grid.
  !> Be careful about doing this. Matrices have pointers to process grids. If
  !> you destruct a process grid without destructing the matrices pointing to
  !> it, they will become unusable.
  RECURSIVE SUBROUTINE DestructProcessGrid(grid_in)
    !> The grid to destruct. If none specified this destroys the global grid.
    TYPE(ProcessGrid_t), OPTIONAL, INTENT(INOUT) :: grid_in
    !! Counters
    INTEGER :: II, JJ
    INTEGER :: Ierr

    !! Handle optional parameters
    IF (.NOT. PRESENT(grid_in)) THEN
       CALL DestructProcessGrid(global_grid)
    ELSE !! Destruct
       IF (ALLOCATED(grid_in%blocked_row_comm)) THEN
          CALL MPI_COMM_FREE(grid_in%global_comm, ierr)
          CALL MPI_COMM_FREE(grid_in%row_comm, ierr)
          CALL MPI_COMM_FREE(grid_in%column_comm, ierr)
          CALL MPI_COMM_FREE(grid_in%within_slice_comm, ierr)
          CALL MPI_COMM_FREE(grid_in%between_slice_comm, ierr)
          DO II = 1, grid_in%number_of_blocks_rows
             CALL MPI_COMM_FREE(grid_in%blocked_row_comm(II), ierr)
          END DO
          DEALLOCATE(grid_in%blocked_row_comm)
       END IF

       IF (ALLOCATED(grid_in%blocked_column_comm)) THEN
          DO JJ = 1, grid_in%number_of_blocks_columns
             CALL MPI_COMM_FREE(grid_in%blocked_column_comm(JJ), ierr)
          END DO
          DEALLOCATE(grid_in%blocked_column_comm)
       END IF

       IF (ALLOCATED(grid_in%blocked_within_slice_comm)) THEN
          DO JJ = 1, grid_in%number_of_blocks_columns
             DO II = 1, grid_in%number_of_blocks_rows
                CALL MPI_COMM_FREE(grid_in%blocked_within_slice_comm(II, JJ), &
                     & ierr)
             END DO
          END DO
          DEALLOCATE(grid_in%blocked_within_slice_comm)
       END IF

       IF (ALLOCATED(grid_in%blocked_between_slice_comm)) THEN
          DO JJ = 1, grid_in%number_of_blocks_columns
             DO II = 1, grid_in%number_of_blocks_rows
                CALL MPI_COMM_FREE(grid_in%blocked_between_slice_comm(II, JJ), &
                     & ierr)
             END DO
          END DO
          DEALLOCATE(grid_in%blocked_between_slice_comm)
       END IF
    END IF

  END SUBROUTINE DestructProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a process grid, this splits it into two grids of even size
  SUBROUTINE SplitProcessGrid(old_grid, new_grid, my_color, split_slice, &
       & between_grid_comm)
    !> The old grid to split
    TYPE(ProcessGrid_t), INTENT(INOUT) :: old_grid
    !> The new grid that we are creating
    TYPE(ProcessGrid_t), INTENT(INOUT) :: new_grid
    !> A color value indicating which set this process was split into
    INTEGER, INTENT(OUT) :: my_color
    !> True if we were able to split along slices.
    LOGICAL, INTENT(OUT) :: split_slice
    !> A communicator for sending messages between groups.
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
       midpoint = old_grid%num_process_slices / 2
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
       left_grid_size = midpoint * cols * rows
       !! Next try to split the bigger direction
    ELSE IF (old_grid%num_process_rows .GT. old_grid%num_process_columns) THEN
       midpoint = old_grid%num_process_rows / 2
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
       left_grid_size = midpoint * cols * slices
       !! Default Case
    ELSE
       midpoint = old_grid%num_process_columns / 2
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
       left_grid_size = midpoint * slices * rows
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
  FUNCTION IsRoot(grid_in) RESULT(is_root)
    !! Parameters
    !> The process grid.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid_in
    !> True if the current process is root.
    LOGICAL :: is_root

    IF (PRESENT(grid_in)) THEN
       is_root = (grid_in%global_rank .EQ. 0)
    ELSE
       is_root = (global_grid%global_rank .EQ. 0)
    END IF
  END FUNCTION IsRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the slice of the current process.
  FUNCTION GetMySlice(grid) RESULT(return_val)
    !> The process grid.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid
    !> Slice number of the current process.
    INTEGER :: return_val

    IF (PRESENT(grid)) THEN
       return_val = grid%my_slice
    ELSE
       return_val = global_grid%my_slice
    END IF
  END FUNCTION GetMySlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the column of the current process.
  FUNCTION GetMyColumn(grid) RESULT(return_val)
    !> The process grid.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid
    !> The column number of the current process.
    INTEGER :: return_val

    IF (PRESENT(grid)) THEN
       return_val = grid%my_column
    ELSE
       return_val = global_grid%my_column
    END IF
  END FUNCTION GetMyColumn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the row of the current process.
  FUNCTION GetMyRow(grid) RESULT(return_val)
    !> The process grid.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: grid
    !> The row number of the current process.
    INTEGER :: return_val

    IF (PRESENT(grid)) THEN
       return_val = grid%my_row
    ELSE
       return_val = global_grid%my_row
    END IF
  END FUNCTION GetMyRow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sometimes we only want to specify for a process grid the number of slices
  !> and then automatically compute the right number of rows and columns.
  SUBROUTINE ComputeGridSize(total_processors, set_slices, rows, columns)
    !> Total processors in the grid
    INTEGER, INTENT(IN) :: total_processors
    !> Desired number of slices
    INTEGER, INTENT(IN) :: set_slices
    !> Computed number of rows
    INTEGER, INTENT(OUT) :: rows
    !> Computed number of columns
    INTEGER, INTENT(OUT) :: columns
    !! Local variables
    INTEGER :: slice_size, size_search
    INTEGER :: II

    rows   = 1
    columns = 1
    slice_size = total_processors / set_slices
    size_search = FLOOR(SQRT(REAL(slice_size)))
    DO II = size_search, 1, -1
       IF (MOD(slice_size, II) .EQ. 0) THEN
          rows = II
          columns = slice_size / II
          EXIT
       END IF
    END DO

  END SUBROUTINE ComputeGridSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pick an appropriate number of process slices for this calculation.
  !> This routine will focus on whether we can make a valid process grid with
  !> several slices.
  SUBROUTINE ComputeNumSlices(total_processors, slices)
    !> Total processors in the grid.
    INTEGER, INTENT(IN) :: total_processors
    !> Number of slices to use.
    INTEGER, INTENT(OUT) :: slices
    !! Local Variables
    INTEGER :: slice_size
    INTEGER :: slice_dim
    LOGICAL :: found

    !! Try manually values [4, 3, 2]. If they do not work, give up and use 1.
    found = .FALSE.
    DO slices = MIN(4, total_processors), 2, -1
       slice_size = total_processors / slices
       IF (slice_size * slices .NE. total_processors) CYCLE

       !! First try a square grid.
       slice_dim  = FLOOR(SQRT(REAL(slice_size)))
       IF (slice_dim*slice_dim .EQ. slice_size) THEN
          FOUND = .TRUE.
          EXIT
       END IF

       !! If not, we try a grid where the rows are twice the number of columns.
       slice_dim = FLOOR(SQRT(REAL(slice_size/2)))
       IF (slice_dim*slice_dim * 2 .EQ. slice_size) THEN
          FOUND = .TRUE.
          EXIT
       END IF
    END DO
    IF (.NOT. FOUND) slices = 1

  END SUBROUTINE ComputeNumSlices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write out some basic information about this process grid to the log.
  RECURSIVE SUBROUTINE WriteProcessGridInfo(this)
    !> The grid to print about. If not specified, global information printed.
    TYPE(ProcessGrid_t), OPTIONAL, INTENT(IN) :: this

    IF (PRESENT(this)) THEN
       CALL WriteHeader("Process Grid")
       CALL EnterSubLog
       CALL WriteListElement(key = "Process Rows", &
            & VALUE = this%num_process_rows)
       CALL WriteListElement(key = "Process Columns", &
            & VALUE = this%num_process_columns)
       CALL WriteListElement(key = "Process Slices", &
            & VALUE = this%num_process_slices)
       CALL WriteListElement(key = "Column Blocks", &
            & VALUE = this%number_of_blocks_columns)
       CALL WriteListElement(key = "Row Blocks", &
            & VALUE = this%number_of_blocks_rows)
       CALL ExitSubLog
    ELSE
       CALL WriteProcessGridInfo(global_grid)
    END IF
  END SUBROUTINE WriteProcessGridInfo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ProcessGridModule
