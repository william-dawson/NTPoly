!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Operations.
MODULE MatrixPSModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE LoggingModule, ONLY : &
       & EnterSubLog, ExitSubLog, WriteElement, WriteListElement, WriteHeader
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & ReduceAndSumMatrixData, ReduceAndSumMatrixCleanup, &
       & TestReduceSizeRequest, TestReduceOuterRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE MatrixMarketModule
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation
  USE ProcessGridModule, ONLY : ProcessGrid_t, global_grid, IsRoot, &
       & SplitProcessGrid, CopyProcessGrid
  USE MatrixSModule
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletModule, ONLY : Triplet_r, Triplet_c, GetMPITripletType_r
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, &
       & DestructTripletList, SortTripletList, AppendToTripletList, &
       & SymmetrizeTripletList, GetTripletAt, RedistributeTripletLists, &
       & ShiftTripletList
  USE iso_c_binding
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for a distributed blocked CSR matrix.
  TYPE, PUBLIC :: Matrix_ps
     !> Number of matrix rows/columns for full matrix, scaled for process grid.
     INTEGER :: logical_matrix_dimension
     !> Number of matrix rows/columns for the full matrix, unscaled.
     INTEGER :: actual_matrix_dimension
     !! Local Storage
     !> A 2D array of local CSR matrices.
     TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: local_data_r
     !> A 2D array of local CSC matrices.
     TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: local_data_c
     INTEGER :: start_column !< first column stored locally.
     INTEGER :: end_column !< last column stored locally  is less than this.
     INTEGER :: start_row !< first row stored locally.
     INTEGER :: end_row !< last row stored locally is less than this.
     INTEGER :: local_columns !< number of local columns.
     INTEGER :: local_rows !< number of local rows.
     TYPE(ProcessGrid_t) :: process_grid !< process grid to operate on
     LOGICAL :: is_complex !< true if the matrix data is true.
  END TYPE Matrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructors/Destructors
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: DestructMatrix
  PUBLIC :: CopyMatrix
  !! File I/O
  PUBLIC :: ConstructMatrixFromMatrixMarket
  PUBLIC :: ConstructMatrixFromBinary
  PUBLIC :: WriteMatrixToMatrixMarket
  PUBLIC :: WriteMatrixToBinary
  !! Fill In Special Matrices
  PUBLIC :: FillMatrixFromTripletList
  PUBLIC :: FillMatrixIdentity
  PUBLIC :: FillMatrixPermutation
  !! Basic Accessors
  PUBLIC :: GetMatrixActualDimension
  PUBLIC :: GetMatrixLogicalDimension
  PUBLIC :: GetMatrixTripletList
  PUBLIC :: GetMatrixBlock
  !! Printing To The Console
  PUBLIC :: PrintMatrix
  PUBLIC :: PrintMatrixInformation
  !! Utilities
  PUBLIC :: GetMatrixLoadBalance
  PUBLIC :: GetMatrixSize
  PUBLIC :: FilterMatrix
  PUBLIC :: MergeMatrixLocalBlocks
  PUBLIC :: SplitMatrixToLocalBlocks
  PUBLIC :: TransposeMatrix
  PUBLIC :: CommSplitMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrix_ps
     MODULE PROCEDURE ConstructEmptyMatrix_ps_cp
  END INTERFACE
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_ps
  END INTERFACE
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_ps
  END INTERFACE
  INTERFACE ConstructMatrixFromMatrixMarket
     MODULE PROCEDURE ConstructMatrixFromMatrixMarket_ps
  END INTERFACE
  INTERFACE ConstructMatrixFromBinary
     MODULE PROCEDURE ConstructMatrixFromBinary_ps
  END INTERFACE
  INTERFACE WriteMatrixToMatrixMarket
     MODULE PROCEDURE WriteMatrixToMatrixMarket_ps
  END INTERFACE
  INTERFACE WriteMatrixToBinary
     MODULE PROCEDURE WriteMatrixToBinary_ps
  END INTERFACE
  INTERFACE FillMatrixFromTripletList
     MODULE PROCEDURE FillMatrixFromTripletList_psr
     MODULE PROCEDURE FillMatrixFromTripletList_psc
  END INTERFACE
  INTERFACE FillMatrixIdentity
     MODULE PROCEDURE FillMatrixIdentity_ps
  END INTERFACE
  INTERFACE FillMatrixPermutation
     MODULE PROCEDURE FillMatrixPermutation_ps
  END INTERFACE
  INTERFACE GetMatrixActualDimension
     MODULE PROCEDURE GetMatrixActualDimension_ps
  END INTERFACE
  INTERFACE GetMatrixLogicalDimension
     MODULE PROCEDURE GetMatrixLogicalDimension_ps
  END INTERFACE
  INTERFACE GetMatrixTripletList
     MODULE PROCEDURE GetMatrixTripletList_psr
     MODULE PROCEDURE GetMatrixTripletList_psc
  END INTERFACE
  INTERFACE GetMatrixBlock
     MODULE PROCEDURE GetMatrixBlock_psr
     MODULE PROCEDURE GetMatrixBlock_psc
  END INTERFACE
  INTERFACE PrintMatrix
     MODULE PROCEDURE PrintMatrix_ps
  END INTERFACE
  INTERFACE PrintMatrixInformation
     MODULE PROCEDURE PrintMatrixInformation_ps
  END INTERFACE
  INTERFACE GetMatrixLoadBalance
     MODULE PROCEDURE GetMatrixLoadBalance_ps
  END INTERFACE
  INTERFACE GetMatrixSize
     MODULE PROCEDURE GetMatrixSize_ps
  END INTERFACE
  INTERFACE FilterMatrix
     MODULE PROCEDURE FilterMatrix_ps
  END INTERFACE
  INTERFACE RedistributeData
     MODULE PROCEDURE RedistributeData_psr
     MODULE PROCEDURE RedistributeData_psc
  END INTERFACE
  INTERFACE MergeMatrixLocalBlocks
     MODULE PROCEDURE MergeMatrixLocalBlocks_psr
     MODULE PROCEDURE MergeMatrixLocalBlocks_psc
  END INTERFACE
  INTERFACE SplitMatrixToLocalBlocks
     MODULE PROCEDURE SplitMatrixToLocalBlocks_psr
     MODULE PROCEDURE SplitMatrixToLocalBlocks_psc
  END INTERFACE
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_ps
  END INTERFACE
  INTERFACE CommSplitMatrix
     MODULE PROCEDURE CommSplitMatrix_ps
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix.
  !! @param[out] this the matrix to be constructed.
  !! @param[in] matrix_dim_ the dimension of the full matrix.
  !! @param[in] process_grid_in a process grid to host the matrix (optional).
  !! @param[in] is_complex_in true if you want to use complex numbers (optional)
  SUBROUTINE ConstructEmptyMatrix_ps(this, matrix_dim_, process_grid_in, &
       & is_complex_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT)            :: this
    INTEGER, INTENT(IN)                       :: matrix_dim_
    LOGICAL, INTENT(IN), OPTIONAL             :: is_complex_in
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: process_grid_in
    !! Local Variables
    TYPE(Matrix_lsr) :: zeromatrix_r
    TYPE(Matrix_lsc) :: zeromatrix_c

    CALL DestructMatrix(this)

    !! Process Grid
    IF (PRESENT(process_grid_in)) THEN
       CALL CopyProcessGrid(process_grid_in, this%process_grid)
    ELSE
       CALL CopyProcessGrid(global_grid, this%process_grid)
    END IF

    !! Complex determination
    IF (PRESENT(is_complex_in)) THEN
       this%is_complex = is_complex_in
    ELSE
       this%is_complex = .FALSE.
    END IF

    !! Matrix Dimensions
    this%actual_matrix_dimension = matrix_dim_
    this%logical_matrix_dimension = CalculateScaledDimension(this, matrix_dim_)

    !! Full Local Data Size Description
    this%local_rows = &
         & this%logical_matrix_dimension/this%process_grid%num_process_rows
    this%local_columns = &
         & this%logical_matrix_dimension/this%process_grid%num_process_columns

    !! Which Block Does This Process Hold?
    this%start_row = this%local_rows * this%process_grid%my_row + 1
    this%end_row   = this%start_row + this%local_rows
    this%start_column = this%local_columns * this%process_grid%my_column + 1
    this%end_column   = this%start_column + this%local_columns

    !! Build local storage
    IF (this%is_complex) THEN
       ALLOCATE(this%local_data_c(this%process_grid%number_of_blocks_rows, &
            & this%process_grid%number_of_blocks_columns))
       zeromatrix_c = Matrix_lsc(this%local_rows, this%local_columns)
       CALL SplitMatrixToLocalBlocks(this, zeromatrix_c)
       CALL DestructMatrix(zeromatrix_c)
    ELSE
       ALLOCATE(this%local_data_r(this%process_grid%number_of_blocks_rows, &
            & this%process_grid%number_of_blocks_columns))
       zeromatrix_r = Matrix_lsr(this%local_rows, this%local_columns)
       CALL SplitMatrixToLocalBlocks(this, zeromatrix_r)
       CALL DestructMatrix(zeromatrix_r)
    END IF
  END SUBROUTINE ConstructEmptyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix using another matrix
  !! to determine the parameters. Note that no data is copied, the matrix
  !! will be empty.
  !! @param[out] this the matrix to be constructed.
  !! @param[in] matrix_dim_ the dimension of the full matrix.
  !! @param[in] process_grid_in a process grid to host the matrix (optional).
  !! @param[in] is_complex_in true if you want to use complex numbers (optional)
  SUBROUTINE ConstructEmptyMatrix_ps_cp(this, reference_matrix)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: reference_matrix

    CALL ConstructEmptyMatrix(this, reference_matrix%actual_matrix_dimension, &
         & reference_matrix%process_grid, reference_matrix%is_complex)
  END SUBROUTINE ConstructEmptyMatrix_ps_cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix.
  !! @param[inout] this the matrix to destruct.
  PURE SUBROUTINE DestructMatrix_ps(this)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: II, JJ

    IF (ALLOCATED(this%local_data_r)) THEN
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL DestructMatrix(this%local_data_r(II,JJ))
             CALL DestructMatrix(this%local_data_c(II,JJ))
          END DO
       END DO
       DEALLOCATE(this%local_data_r)
       DEALLOCATE(this%local_data_c)
    END IF
  END SUBROUTINE DestructMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopyMatrix_ps(matA, matB)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(INOUT) :: matB

    CALL DestructMatrix(matB)
    matB = matA
  END SUBROUTINE CopyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  !! Read \cite boisvert1996matrix for the details.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructMatrixFromMatrixMarket_ps(this, file_name)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    INTEGER, PARAMETER :: MAX_LINE_LENGTH = 100
    !! File Handles
    INTEGER :: local_file_handler
    INTEGER :: mpi_file_handler
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Reading The File
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: matrix_rows, matrix_columns, total_values
    !! Length Variables
    INTEGER :: header_length
    INTEGER(KIND=MPI_OFFSET_KIND) :: total_file_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_data_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_data_size_plus_buffer
    INTEGER :: current_line_length
    !! Input Buffers
    CHARACTER(len=MAX_LINE_LENGTH) :: input_buffer
    CHARACTER(len=:), ALLOCATABLE :: mpi_input_buffer
    CHARACTER(len=MAX_LINE_LENGTH) :: temp_substring
    !! Temporary Variables
    INTEGER :: bytes_per_character
    CHARACTER(len=1) :: temp_char
    LOGICAL :: found_comment_line
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER :: full_buffer_counter
    LOGICAL :: end_of_buffer
    LOGICAL :: error_occured
    INTEGER :: ierr

    !! Setup Involves Just The Root Opening And Reading Parameter Data
    CALL StartTimer("MPI Read Text")
    bytes_per_character = sizeof(temp_char)
    IF (IsRoot(global_grid)) THEN
       header_length = 0
       local_file_handler = 16
       OPEN(local_file_handler, file=file_name, iostat=ierr, status="old")
       IF (ierr .EQ. 0) THEN
          !! Parse the header.
          READ(local_file_handler,fmt='(A)') input_buffer
          error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
               & pattern_type)
          header_length = header_length + LEN_TRIM(input_buffer) + 1
          !! First Read In The Comment Lines
          found_comment_line = .TRUE.
          DO WHILE (found_comment_line)
             READ(local_file_handler,fmt='(A)') input_buffer
             !! +1 for newline
             header_length = header_length + LEN_TRIM(input_buffer) + 1
             IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
                found_comment_line = .FALSE.
             END IF
          END DO
          !! Get The Matrix Parameters
          READ(input_buffer,*) matrix_rows, matrix_columns, total_values
          CLOSE(local_file_handler)
       ELSE
          WRITE(*,*) file_name, " doesn't exist"
       END IF
    ELSE
       ierr = 0
    END IF

    IF (ierr .NE. 0) THEN
       CALL MPI_Abort(global_grid%global_comm, -1, ierr)
    END IF

    !! Broadcast Parameters
    CALL MPI_Bcast(matrix_rows, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(matrix_columns, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(total_values, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(header_length, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(sparsity_type, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(data_type, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(pattern_type, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)

    !! Build Local Storage
    CALL ConstructEmptyMatrix(this, matrix_rows, global_grid)

    !! Global read
    CALL MPI_File_open(this%process_grid%global_comm,file_name,MPI_MODE_RDONLY,&
         & MPI_INFO_NULL,mpi_file_handler,ierr)
    CALL MPI_File_get_size(mpi_file_handler,total_file_size,ierr)

    !! Compute Offsets And Data Size
    local_data_size = (total_file_size - bytes_per_character*header_length)/&
         & this%process_grid%total_processors
    IF (local_data_size .LT. 2*MAX_LINE_LENGTH) THEN
       local_data_size = 2*MAX_LINE_LENGTH
    END IF
    local_offset = bytes_per_character*header_length + &
         local_data_size*this%process_grid%global_rank

    !! Check if this processor has any work to do, and set the appropriate
    !! buffer size. We also add some buffer space, so you can read beyond
    !! your local data size in case the local data read ends in the middle
    !! of a line.
    IF (local_offset .LT. total_file_size) THEN
       local_data_size_plus_buffer = local_data_size + &
            & MAX_LINE_LENGTH*bytes_per_character
       IF (local_offset + local_data_size_plus_buffer .GT. total_file_size) THEN
          local_data_size_plus_buffer = (total_file_size - local_offset)
       END IF
    ELSE
       local_data_size_plus_buffer = 0
    END IF

    !! A buffer to read the data into.
    ALLOCATE(CHARACTER(LEN=local_data_size_plus_buffer) :: mpi_input_buffer)

    !! Do Actual Reading
    CALL MPI_File_read_at_all(mpi_file_handler,local_offset,mpi_input_buffer, &
         & INT(local_data_size_plus_buffer),MPI_CHARACTER,mpi_status,ierr)

    !! Trim Off The Half Read Line At The Start
    IF (.NOT. this%process_grid%global_rank .EQ. this%process_grid%RootID) THEN
       full_buffer_counter = INDEX(mpi_input_buffer,new_line('A')) + 1
    ELSE
       full_buffer_counter = 1
    END IF

    !! Read By Line
    end_of_buffer = .FALSE.
    IF (local_data_size_plus_buffer .EQ. 0) THEN
       end_of_buffer = .TRUE.
    END IF
    triplet_list = TripletList_r()
    DO WHILE(.NOT. end_of_buffer)
       current_line_length = INDEX(mpi_input_buffer(full_buffer_counter:),&
            new_line('A'))

       IF (current_line_length .EQ. 0) THEN !! Hit The End Of The Buffer
          end_of_buffer = .TRUE.
       ELSE
          temp_substring = mpi_input_buffer( &
               & full_buffer_counter:full_buffer_counter+current_line_length-1)
          IF (current_line_length .GT. 1) THEN
             READ(temp_substring(:current_line_length-1),*) &
                  & temp_triplet%index_row, temp_triplet%index_column, &
                  & temp_triplet%point_value
             CALL AppendToTripletList(triplet_list, temp_triplet)
          END IF

          IF (full_buffer_counter + current_line_length .GE. &
               & local_data_size+2) THEN
             IF (.NOT. this%process_grid%global_rank .EQ. &
                  & this%process_grid%total_processors-1) THEN
                end_of_buffer = .TRUE.
             END IF
          END IF
          full_buffer_counter = full_buffer_counter + current_line_length
       END IF
    END DO

    !! Cleanup
    CALL MPI_File_close(mpi_file_handler,ierr)
    CALL StopTimer("MPI Read Text")
    CALL MPI_Barrier(this%process_grid%global_comm,ierr)

    !! Redistribute The Matrix
    CALL SymmetrizeTripletList(triplet_list, pattern_type)
    CALL FillMatrixFromTripletList(this,triplet_list)

    CALL DestructTripletList(triplet_list)
    DEALLOCATE(mpi_input_buffer)
  END SUBROUTINE ConstructMatrixFromMatrixMarket_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  !! Faster than text, so this is good for check pointing.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructMatrixFromBinary_ps(this, file_name)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! File Handles
    INTEGER :: mpi_file_handler
    !! Reading The File
    TYPE(TripletList_r) :: triplet_list
    INTEGER :: matrix_rows, matrix_columns, total_values
    INTEGER, DIMENSION(3) :: matrix_information
    INTEGER :: local_triplets
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER :: bytes_per_int, bytes_per_double
    INTEGER :: triplet_mpi_type
    !! Temporary variables
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER :: ierr

    CALL StartTimer("MPI Read Binary")
    CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
    CALL MPI_Type_extent(MPINTREAL,bytes_per_double,ierr)

    CALL MPI_File_open(global_grid%global_comm,file_name,MPI_MODE_RDONLY,&
         & MPI_INFO_NULL,mpi_file_handler,ierr)
    IF (ierr .NE. 0) THEN
       IF (IsRoot(global_grid)) THEN
          WRITE(*,*) file_name, " doesn't exist"
       END IF
       CALL MPI_Abort(global_grid%global_comm, -1, ierr)
    END IF

    !! Get The Matrix Parameters
    IF (IsRoot(global_grid)) THEN
       local_offset = 0
       CALL MPI_File_read_at(mpi_file_handler, local_offset, &
            & matrix_information, 3, MPI_INT, mpi_status, ierr)
       matrix_rows = matrix_information(1)
       matrix_columns = matrix_information(2)
       total_values = matrix_information(3)
    END IF

    !! Broadcast Parameters
    CALL MPI_Bcast(matrix_rows, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(matrix_columns, 1, MPI_INT, global_grid%RootID, &
         & global_grid%global_comm, ierr)
    CALL MPI_Bcast(total_values, 1, MPI_INT ,global_grid%RootID, &
         & global_grid%global_comm, ierr)

    !! Build Local Storage
    CALL ConstructEmptyMatrix(this, matrix_rows, global_grid)

    !! Compute Offset
    local_triplets = total_values/this%process_grid%total_processors
    local_offset = local_triplets * (this%process_grid%global_rank)
    header_size = 3 * bytes_per_int
    IF (this%process_grid%global_rank .EQ. &
         & this%process_grid%total_processors - 1) THEN
       local_triplets = INT(total_values) - INT(local_offset)
    END IF
    local_offset = local_offset*(bytes_per_int*2+bytes_per_double) + header_size
    triplet_list = TripletList_r(local_triplets)

    !! Do The Actual Reading
    triplet_mpi_type = GetMPITripletType_r()
    CALL MPI_File_set_view(mpi_file_handler,local_offset,triplet_mpi_type,&
         & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
    CALL MPI_File_read_all(mpi_file_handler,triplet_list%data,local_triplets,&
         & triplet_mpi_type,mpi_status,ierr)
    CALL MPI_File_close(mpi_file_handler,ierr)
    CALL StopTimer("MPI Read Binary")

    CALL FillMatrixFromTripletList(this,triplet_list)

    !! Cleanup
    CALL DestructTripletList(triplet_list)
  END SUBROUTINE ConstructMatrixFromBinary_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a binary file.
  !! Faster than text, so this is good for check pointing.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteMatrixToBinary_ps(this,file_name)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
    INTEGER :: mpi_file_handler
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_data_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: write_offset
    !! For The Special Datatype
    INTEGER :: triplet_mpi_type
    !! Temporary Variables
    INTEGER :: temp_int
    REAL(NTREAL) :: temp_double
    INTEGER :: bytes_per_int, bytes_per_double
    INTEGER, DIMENSION(3) :: header_buffer
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: zero_offset = 0
    INTEGER :: counter
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(this, merged_local_data)

    !! Determine Write Location
    bytes_per_int = sizeof(temp_int)
    bytes_per_double = sizeof(temp_double)
    local_data_size = SIZE(merged_local_data%values)*(bytes_per_int*2 + &
         & bytes_per_double*1)
    header_size = bytes_per_int*3
    ALLOCATE(local_values_buffer(this%process_grid%slice_size))
    CALL MPI_Allgather(SIZE(merged_local_data%values),1,MPI_INT,&
         & local_values_buffer,1,MPI_INT,&
         & this%process_grid%within_slice_comm,ierr)
    write_offset = 0
    write_offset = write_offset + header_size
    DO counter = 1,this%process_grid%within_slice_rank
       write_offset = write_offset + &
            & local_values_buffer(counter)*(bytes_per_int*2+bytes_per_double*1)
    END DO

    !! Write The File
    IF (this%process_grid%between_slice_rank .EQ. 0) THEN
       !! Create Special MPI Type
       triplet_mpi_type = GetMPITripletType_r()

       CALL MatrixToTripletList(merged_local_data, triplet_list)
       !! Absolute Positions
       CALL ShiftTripletList(triplet_list, this%start_row - 1, &
            & this%start_column - 1)
       CALL MPI_File_open(this%process_grid%within_slice_comm,file_name,&
            & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL, &
            & mpi_file_handler, ierr)
       !! Write Header
       IF (this%process_grid%within_slice_rank .EQ. 0) THEN
          header_buffer(1) = this%actual_matrix_dimension
          header_buffer(2) = this%actual_matrix_dimension
          header_buffer(3) = SUM(local_values_buffer)
          CALL MPI_File_write_at(mpi_file_handler,zero_offset,header_buffer,3,&
               & MPI_INT,mpi_status,ierr)
       END IF
       !! Write The Rest
       CALL MPI_File_set_view(mpi_file_handler,write_offset,triplet_mpi_type,&
            & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
       CALL MPI_File_write(mpi_file_handler,triplet_list%data, &
            & triplet_list%CurrentSize, triplet_mpi_type,MPI_STATUS_IGNORE, &
            & ierr)

       !! Cleanup
       CALL MPI_File_close(mpi_file_handler,ierr)
       CALL DestructTripletList(triplet_list)
    END IF
    DEALLOCATE(local_values_buffer)
    CALL MPI_Barrier(this%process_grid%global_comm,ierr)
    CALL DestructMatrix(merged_local_data)
  END SUBROUTINE WriteMatrixToBinary_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a distributed sparse matrix to a matrix market file.
  !! Read \cite boisvert1996matrix for the details.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteMatrixToMatrixMarket_ps(this,file_name)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    INTEGER, PARAMETER :: MAX_LINE_LENGTH = 1024
    !! Local MPI Variables
    INTEGER :: mpi_file_handler
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    INTEGER :: triplet_list_string_length
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: write_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_offset
    INTEGER(KIND=MPI_OFFSET_KIND), PARAMETER :: zero_size = 0
    !! Strings
    CHARACTER(len=:), ALLOCATABLE :: header_line1
    CHARACTER(len=:), ALLOCATABLE :: header_line2
    CHARACTER(len=:), ALLOCATABLE :: write_buffer
    !! Temporary Values
    INTEGER :: counter
    INTEGER :: offset_counter
    INTEGER :: NEW_LINE_LENGTH
    CHARACTER(len=MAX_LINE_LENGTH*2) :: temp_string1
    CHARACTER(len=MAX_LINE_LENGTH) :: temp_string2
    INTEGER :: temp_length
    INTEGER :: bytes_per_character
    CHARACTER(len=1) :: temp_char
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(this, merged_local_data)

    bytes_per_character = sizeof(temp_char)

    !! Create the matrix size line
    NEW_LINE_LENGTH = LEN(new_line('A'))
    WRITE(temp_string1,'(A)') "%%MatrixMarket matrix coordinate real general" &
         & //new_line('A')//"%"//new_line('A')
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string1)) :: header_line1)
    header_line1 = TRIM(temp_string1)

    WRITE(temp_string2,*) this%actual_matrix_dimension, &
         & this%actual_matrix_dimension, GetMatrixSize(this)
    !! I don't understand why the +1 is needed, but it is.
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string2)+NEW_LINE_LENGTH+1) :: header_line2)
    WRITE(header_line2,*) TRIM(temp_string2)//new_line('A')

    header_size = LEN(header_line1) + LEN(header_line2)

    !! Local Data
    CALL MatrixToTripletList(merged_local_data, triplet_list)

    !! Absolute Positions
    CALL ShiftTripletList(triplet_list, this%start_row - 1, &
         & this%start_column - 1)

    !! Figure out the length of the string for storing.
    triplet_list_string_length = 0
    DO counter = 1, triplet_list%CurrentSize
       WRITE(temp_string2,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value, &
            & new_line('A')
       triplet_list_string_length = triplet_list_string_length + &
            & LEN_TRIM(temp_string2)
       triplet_list_string_length = triplet_list_string_length + NEW_LINE_LENGTH
    END DO

    !! Write that string to the write buffer
    ALLOCATE(CHARACTER(len=triplet_list_string_length+1) :: write_buffer)
    offset_counter = 1
    DO counter = 1, triplet_list%CurrentSize
       WRITE(temp_string2,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value, &
            & new_line('A')
       temp_length = LEN_TRIM(temp_string2)+NEW_LINE_LENGTH
       WRITE(write_buffer(offset_counter:offset_counter+temp_length),*) &
            & temp_string2(1:temp_length)
       offset_counter = offset_counter + temp_length
    END DO

    !! Figure out the offset sizes
    ALLOCATE(local_values_buffer(this%process_grid%slice_size))
    CALL MPI_Allgather(triplet_list_string_length,1,MPI_INT,&
         & local_values_buffer,1,MPI_INT,this%process_grid%within_slice_comm,&
         & ierr)
    write_offset = 0
    write_offset = write_offset + header_size
    DO counter = 1,this%process_grid%within_slice_rank
       write_offset = write_offset + &
            & local_values_buffer(counter)
    END DO

    !! Global Write
    IF (this%process_grid%between_slice_rank .EQ. 0) THEN
       CALL MPI_File_open(this%process_grid%within_slice_comm,file_name, &
            & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL, &
            & mpi_file_handler,ierr)
       CALL MPI_File_set_size(mpi_file_handler,zero_size,ierr)
       !! Write Header
       IF (this%process_grid%within_slice_rank .EQ. 0) THEN
          header_offset = 0
          CALL MPI_File_write_at(mpi_file_handler,header_offset,header_line1, &
               & LEN(header_line1), MPI_CHARACTER,mpi_status,ierr)
          header_offset = header_offset + LEN(header_line1)
          CALL MPI_File_write_at(mpi_file_handler,header_offset,header_line2, &
               & LEN(header_line2), MPI_CHARACTER,mpi_status,ierr)
       END IF
       !! Write Local Data
       CALL MPI_File_set_view(mpi_file_handler,write_offset,MPI_CHARACTER,&
            & MPI_CHARACTER,"native",MPI_INFO_NULL,ierr)
       CALL MPI_File_write(mpi_file_handler,write_buffer, &
            & triplet_list_string_length,&
            & MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

       !! Cleanup
       CALL MPI_File_close(mpi_file_handler,ierr)
    END IF
    CALL MPI_Barrier(this%process_grid%global_comm,ierr)
  END SUBROUTINE WriteMatrixToMatrixMarket_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !! should pass in triplet lists with global coordinates. It doesn't matter
  !! where each triplet is stored, as long as global coordinates are given.
  !! @param[inout] this the matrix to fill.
  !! @param[in] triplet_list the triplet list of values.
  !! @param[in] preduplicated_in if lists are preduplicated across
  !! slices set this to true (optional, default=False).
  SUBROUTINE FillMatrixFromTripletList_psr(this,triplet_list,preduplicated_in)
    !! Parameters
    TYPE(Matrix_ps) :: this
    TYPE(TripletList_r) :: triplet_list
    LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    !! Local Data
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix
    TYPE(Matrix_lsr) :: gathered_matrix

    INCLUDE "includes/FillMatrixFromTripletList.f90"
  END SUBROUTINE FillMatrixFromTripletList_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !! should pass in triplet lists with global coordinates. It doesn't matter
  !! where each triplet is stored, as long as global coordinates are given.
  !! @param[inout] this the matrix to fill.
  !! @param[in] triplet_list the triplet list of values.
  !! @param[in] preduplicated_in if lists are preduplicated across
  !! slices set this to true (optional, default=False).
  SUBROUTINE FillMatrixFromTripletList_psc(this,triplet_list,preduplicated_in)
    !! Parameters
    TYPE(Matrix_ps) :: this
    TYPE(TripletList_c) :: triplet_list
    LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    !! Local Data
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix
    TYPE(Matrix_lsc) :: gathered_matrix

    INCLUDE "includes/FillMatrixFromTripletList.f90"
  END SUBROUTINE FillMatrixFromTripletList_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] this the matrix being filled.
  SUBROUTINE FillMatrixIdentity_ps(this)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this

    IF (this%is_complex) THEN
       CALL FillMatrixIdentity_psc(this)
    ELSE
       CALL FillMatrixIdentity_psr(this)
    END IF

  END SUBROUTINE FillMatrixIdentity_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] this the matrix being filled.
  SUBROUTINE FillMatrixIdentity_psr(this)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix

    INCLUDE "includes/FillMatrixIdentity.f90"

  END SUBROUTINE FillMatrixIdentity_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] this the matrix being filled.
  SUBROUTINE FillMatrixIdentity_psc(this)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: unsorted_triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix

    INCLUDE "includes/FillMatrixIdentity.f90"

  END SUBROUTINE FillMatrixIdentity_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !! If you don't specify permuterows, will default to permuting rows.
  !! @param[inout] this the matrix being filled.
  !! @param[in] permutation_vector describes for each row/column, where it goes.
  !! @param[in] permute_rows_in if true permute rows, false permute columns.
  SUBROUTINE FillMatrixPermutation_ps(this, permutation_vector, permute_rows_in)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    LOGICAL, OPTIONAL, INTENT(IN) :: permute_rows_in
    !! Local Data
    LOGICAL :: permute_rows

    !! Figure out what type of permutation
    IF (PRESENT(permute_rows_in) .AND. permute_rows_in .EQV. .FALSE.) THEN
       permute_rows = .FALSE.
    ELSE
       permute_rows = .TRUE.
    END IF

    IF (this%is_complex) THEN
       CALL FillMatrixPermutation_psc(this, permutation_vector, permute_rows)
    ELSE
       CALL FillMatrixPermutation_psr(this, permutation_vector, permute_rows)
    END IF

  END SUBROUTINE FillMatrixPermutation_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !! If you don't specify permuterows, will default to permuting rows.
  !! @param[inout] this the matrix being filled.
  !! @param[in] permutation_vector describes for each row/column, where it goes.
  !! @param[in] permute_rows_in if true permute rows, false permute columns.
  SUBROUTINE FillMatrixPermutation_psr(this, permutation_vector, rows)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    LOGICAL, INTENT(IN) :: rows
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix

    INCLUDE "includes/FillMatrixPermutation.f90"

  END SUBROUTINE FillMatrixPermutation_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !! If you don't specify permuterows, will default to permuting rows.
  !! @param[inout] this the matrix being filled.
  !! @param[in] permutation_vector describes for each row/column, where it goes.
  !! @param[in] permute_rows_in if true permute rows, false permute columns.
  SUBROUTINE FillMatrixPermutation_psc(this, permutation_vector, rows)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    LOGICAL, INTENT(IN) :: rows
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: unsorted_triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix

    INCLUDE "includes/FillMatrixPermutation.f90"

  END SUBROUTINE FillMatrixPermutation_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  PURE SUBROUTINE GetMatrixTripletList_psr(this, triplet_list)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !! Local Data
    TYPE(Matrix_lsr) :: merged_local_data

    INCLUDE "includes/GetMatrixTripletList.f90"
  END SUBROUTINE GetMatrixTripletList_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  PURE SUBROUTINE GetMatrixTripletList_psc(this, triplet_list)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !! Local Data
    TYPE(Matrix_lsc) :: merged_local_data

    INCLUDE "includes/GetMatrixTripletList.f90"
  END SUBROUTINE GetMatrixTripletList_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list. Block is
  !! defined by the row/column start/end values.
  !! This is slower than GetMatrixTripletList, because communication is required.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  !! @param[in] start_row the starting row for data to store on this process.
  !! @param[in] end_row the ending row for data to store on this process.
  !! @param[in] start_column the starting col for data to store on this process
  !! @param[in] end_column the ending col for data to store on this process
  SUBROUTINE GetMatrixBlock_psr(this, triplet_list, start_row, end_row, &
       & start_column, end_column)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    INTEGER :: start_row, end_row
    INTEGER :: start_column, end_column
    !! Local Data
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(TripletList_r) :: local_triplet_list
    !! Send Buffer
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    !! Receive Buffer
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! Temp Values
    TYPE(Triplet_r) :: temp_triplet

#define MPIDATATYPE MPINTREAL
#include "includes/GetMatrixBlock.f90"
#undef MPIDATATYPE

  END SUBROUTINE GetMatrixBlock_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list. Block is
  !! defined by the row/column start/end values.
  !! This is slower than GetMatrixTripletList, because communication is required.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  !! @param[in] start_row the starting row for data to store on this process.
  !! @param[in] end_row the ending row for data to store on this process.
  !! @param[in] start_column the starting col for data to store on this process
  !! @param[in] end_column the ending col for data to store on this process
  SUBROUTINE GetMatrixBlock_psc(this, triplet_list, start_row, end_row, &
       & start_column, end_column)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    INTEGER :: start_row, end_row
    INTEGER :: start_column, end_column
    !! Local Data
    TYPE(Matrix_lsc) :: merged_local_data
    TYPE(TripletList_c) :: local_triplet_list
    !! Send Buffer
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    !! Receive Buffer
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! Temp Values
    TYPE(Triplet_c) :: temp_triplet

#define MPIDATATYPE MPINTCOMPLEX
#include "includes/GetMatrixBlock.f90"
#undef MPIDATATYPE

  END SUBROUTINE GetMatrixBlock_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the actual dimension of the matrix.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetMatrixActualDimension_ps(this) RESULT(DIMENSION)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%actual_matrix_dimension
  END FUNCTION GetMatrixActualDimension_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the logical dimension of the matrix.
  !! Includes padding.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetMatrixLogicalDimension_ps(this) RESULT(DIMENSION)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%logical_matrix_dimension
  END FUNCTION GetMatrixLogicalDimension_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out information about a distributed sparse matrix.
  !! Sparsity, and load balancing information.
  !! @param[in] this the matrix to print information about.
  SUBROUTINE PrintMatrixInformation_ps(this)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    CALL GetMatrixLoadBalance(this,min_size,max_size)
    sparsity = REAL(GetMatrixSize(this),KIND=NTREAL) / &
         & (REAL(this%actual_matrix_dimension,KIND=NTREAL)**2)

    CALL WriteHeader("Load_Balance")
    CALL EnterSubLog
    CALL WriteListElement(key="min_size", int_value_in=min_size)
    CALL WriteListElement(key="max_size", int_value_in=max_size)
    CALL ExitSubLog
    CALL WriteElement(key="Dimension",int_value_in=this%actual_matrix_dimension)
    CALL WriteElement(key="Sparsity", float_value_in=sparsity)
  END SUBROUTINE PrintMatrixInformation_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a distributed sparse matrix.
  !! This is a serial print routine, and should probably only be used for debug
  !! purposes.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_ps(this, file_name_in)
    !! Parameters
    TYPE(Matrix_ps) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in

    IF (this%is_complex) THEN
       IF (PRESENT(file_name_in)) THEN
          CALL PrintMatrix_psc(this, file_name_in)
       ELSE
          CALL PrintMatrix_psc(this)
       END IF
    ELSE
       IF (PRESENT(file_name_in)) THEN
          CALL PrintMatrix_psr(this, file_name_in)
       ELSE
          CALL PrintMatrix_psr(this)
       END IF
    END IF
  END SUBROUTINE PrintMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a distributed sparse matrix.
  !! This is a serial print routine, and should probably only be used for debug
  !! purposes.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_psr(this, file_name_in)
    !! Parameters
    TYPE(Matrix_ps) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Temporary Variables
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(Matrix_lsr) :: merged_local_dataT
    TYPE(Matrix_lsr) :: merged_columns
    TYPE(Matrix_lsr) :: merged_columnsT
    TYPE(Matrix_lsr) :: full_gathered

    INCLUDE "includes/PrintMatrix.f90"
  END SUBROUTINE PrintMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a distributed sparse matrix.
  !! This is a serial print routine, and should probably only be used for debug
  !! purposes.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_psc(this, file_name_in)
    !! Parameters
    TYPE(Matrix_ps) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Temporary Variables
    TYPE(Matrix_lsc) :: merged_local_data
    TYPE(Matrix_lsc) :: merged_local_dataT
    TYPE(Matrix_lsc) :: merged_columns
    TYPE(Matrix_lsc) :: merged_columnsT
    TYPE(Matrix_lsc) :: full_gathered

    INCLUDE "includes/PrintMatrix.f90"
  END SUBROUTINE PrintMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !! All (absolute) values below the threshold are set to zero.
  !! @param[inout] this matrix to filter
  !! @param[in] threshold (absolute) values below this are filtered
  SUBROUTINE FilterMatrix_ps(this, threshold)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: threshold

    IF (this%is_complex) THEN
       CALL Filter_Matrix_psc(this, threshold)
    ELSE
       CALL Filter_Matrix_psr(this, threshold)
    END IF
  END SUBROUTINE FilterMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !! All (absolute) values below the threshold are set to zero.
  !! @param[inout] this matrix to filter
  !! @param[in] threshold (absolute) values below this are filtered
  SUBROUTINE FilterMatrix_psr(this, threshold)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary

    INCLUDE "includes/FilterMatrix.f90"
  END SUBROUTINE FilterMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !! All (absolute) values below the threshold are set to zero.
  !! @param[inout] this matrix to filter
  !! @param[in] threshold (absolute) values below this are filtered
  SUBROUTINE FilterMatrix_psc(this, threshold)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: new_list
    TYPE(Triplet_c) :: temporary

    INCLUDE "includes/FilterMatrix.f90"
  END SUBROUTINE FilterMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the total number of non zero entries in the distributed sparse matrix.
  !! @param[in] this the distributed sparse matrix to calculate the non-zero
  !! entries of.
  !! @return the number of non-zero entries in the matrix.
  FUNCTION GetMatrixSize_ps(this) RESULT(total_size)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER(c_long) :: total_size
    !! Local Data
    !integer :: local_size
    REAL(NTREAL) :: local_size
    REAL(NTREAL) :: temp_size
    TYPE(Matrix_lsc) :: merged_local_data_c
    TYPE(Matrix_lsr) :: merged_local_data_r
    INTEGER :: ierr

    !! Merge all the local data
    IF (this%is_complex) THEN
       CALL MergeMatrixLocalBlocks(this, merged_local_data_c)
       local_size = SIZE(merged_local_data_c%values)
       CALL DestructMatrix(merged_local_data_c)
    ELSE
       CALL MergeMatrixLocalBlocks(this, merged_local_data_r)
       local_size = SIZE(merged_local_data_r%values)
       CALL DestructMatrix(merged_local_data_r)
    END IF

    !! Global Sum
    CALL MPI_Allreduce(local_size,temp_size,1,MPINTREAL,MPI_SUM,&
         & this%process_grid%within_slice_comm, ierr)
  END FUNCTION GetMatrixSize_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get a measure of how load balanced this matrix is. For each process, the
  !! number of non-zero entries is calculated. Then, this function returns
  !! the max and min of those values.
  !! @param[in] this The matrix to compute the measure on.
  !! @param[out] min_size the minimum entries contained on a single process.
  !! @param[out] max_size the maximum entries contained on a single process.
  SUBROUTINE GetMatrixLoadBalance_ps(this, min_size, max_size)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: min_size
    INTEGER, INTENT(OUT) :: max_size
    !! Local Data
    INTEGER :: local_size
    TYPE(Matrix_lsc) :: merged_local_data_c
    TYPE(Matrix_lsr) :: merged_local_data_r
    INTEGER :: ierr

    !! Merge all the local data
    IF (this%is_complex) THEN
       CALL MergeMatrixLocalBlocks(this, merged_local_data_c)
       local_size = SIZE(merged_local_data_c%values)
       CALL DestructMatrix(merged_local_data_c)
    ELSE
       CALL MergeMatrixLocalBlocks(this, merged_local_data_r)
       local_size = SIZE(merged_local_data_r%values)
       CALL DestructMatrix(merged_local_data_r)
    END IF

    !! Global Reduce
    CALL MPI_Allreduce(local_size,max_size,1,MPI_INT,MPI_MAX,&
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allreduce(local_size,min_size,1,MPI_INT,MPI_MIN,&
         & this%process_grid%within_slice_comm, ierr)

  END SUBROUTINE GetMatrixLoadBalance_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  !! @param[in] AMat The matrix to transpose.
  !! @param[out] TransMat A^T
  SUBROUTINE TransposeMatrix_ps(AMat, TransMat)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat

    IF (AMat%is_complex) THEN
       CALL TransposeMatrix_psr(AMat, TransMat)
    ELSE
       CALL TransposeMatrix_psr(AMat, TransMat)
    END IF

  END SUBROUTINE TransposeMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  !! @param[in] AMat The matrix to transpose.
  !! @param[out] TransMat A^T
  SUBROUTINE TransposeMatrix_psr(AMat, TransMat)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary, temporary_t

    INCLUDE "includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  !! @param[in] AMat The matrix to transpose.
  !! @param[out] TransMat A^T
  SUBROUTINE TransposeMatrix_psc(AMat, TransMat)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat
    !! Local Variables
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: new_list
    TYPE(Triplet_c) :: temporary, temporary_t

    INCLUDE "includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split the current communicator, and give each group a complete copy of this
  !! @param[in] this the matrix to split.
  !! @param[out] split_mat a copy of the matrix hosted on a small process grid.
  !! @param[out] my_color distinguishes between the two groups (optional).
  !! @param[out] split_slice if we split along the slice direction, this is True
  SUBROUTINE CommSplitMatrix_ps(this, split_mat, my_color, split_slice)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    INTEGER, INTENT(OUT) :: my_color
    LOGICAL, INTENT(OUT) :: split_slice

    IF (this%is_complex) THEN
       CALL CommSplitMatrix_psc(this, split_mat, my_color, split_slice)
    ELSE
       CALL CommSplitMatrix_psr(this, split_mat, my_color, split_slice)
    END IF

  END SUBROUTINE CommSplitMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split the current communicator, and give each group a complete copy of this
  !! @param[in] this the matrix to split.
  !! @param[out] split_mat a copy of the matrix hosted on a small process grid.
  !! @param[out] my_color distinguishes between the two groups (optional).
  !! @param[out] split_slice if we split along the slice direction, this is True
  SUBROUTINE CommSplitMatrix_psr(this, split_mat, my_color, split_slice)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    INTEGER, INTENT(OUT) :: my_color
    LOGICAL, INTENT(OUT) :: split_slice
    !! For Data Redistribution
    TYPE(TripletList_r) :: full_list, new_list
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list

    INCLUDE "includes/CommSplitMatrix.f90"

  END SUBROUTINE CommSplitMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split the current communicator, and give each group a complete copy of this
  !! @param[in] this the matrix to split.
  !! @param[out] split_mat a copy of the matrix hosted on a small process grid.
  !! @param[out] my_color distinguishes between the two groups (optional).
  !! @param[out] split_slice if we split along the slice direction, this is True
  SUBROUTINE CommSplitMatrix_psc(this, split_mat, my_color, split_slice)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    INTEGER, INTENT(OUT) :: my_color
    LOGICAL, INTENT(OUT) :: split_slice
    !! For Data Redistribution
    TYPE(TripletList_c) :: full_list, new_list
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: send_list

    INCLUDE "includes/CommSplitMatrix.f90"

  END SUBROUTINE CommSplitMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute the data in a matrix based on row, column list
  !! This will redistribute the data so that the local data are entries in
  !! the rows and columns list. The order of the row list and column list matter
  !! because local data is filled in the same order.
  !! @param[inout] this the matrix to redistribute
  !! @param[in] index_lookup, describing how data is distributed.
  !! @param[in] reverse_lookup, describing how data is distributed.
  !! @param[in] initial_triplet_list is the current triplet list of global
  !! coordinates
  !! @param[out] sorted_triplet_list returns an allocated triplet list with
  !! local coordinates in sorted order.
  SUBROUTINE RedistributeData_psr(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: index_lookup
    INTEGER, DIMENSION(:), INTENT(IN) :: reverse_index_lookup
    TYPE(TripletList_r), INTENT(IN) :: initial_triplet_list
    TYPE(TripletList_r), INTENT(OUT) :: sorted_triplet_list
    !! Local Data
    TYPE(TripletList_r) :: gathered_list
    TYPE(TripletList_r), DIMENSION(this%process_grid%slice_size) :: &
         & send_triplet_lists
    TYPE(Triplet_r) :: temp_triplet

    INCLUDE "includes/RedistributeData.f90"

  END SUBROUTINE RedistributeData_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute the data in a matrix based on row, column list
  !! This will redistribute the data so that the local data are entries in
  !! the rows and columns list. The order of the row list and column list matter
  !! because local data is filled in the same order.
  !! @param[inout] this the matrix to redistribute
  !! @param[in] index_lookup, describing how data is distributed.
  !! @param[in] reverse_lookup, describing how data is distributed.
  !! @param[in] initial_triplet_list is the current triplet list of global
  !! coordinates
  !! @param[out] sorted_triplet_list returns an allocated triplet list with
  !! local coordinates in sorted order.
  SUBROUTINE RedistributeData_psc(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: index_lookup
    INTEGER, DIMENSION(:), INTENT(IN) :: reverse_index_lookup
    TYPE(TripletList_c), INTENT(IN) :: initial_triplet_list
    TYPE(TripletList_c), INTENT(OUT) :: sorted_triplet_list
    !! Local Data
    TYPE(TripletList_c) :: gathered_list
    TYPE(TripletList_c), DIMENSION(this%process_grid%slice_size) :: &
         & send_triplet_lists
    TYPE(Triplet_c) :: temp_triplet

    INCLUDE "includes/RedistributeData.f90"

  END SUBROUTINE RedistributeData_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate a matrix size that can be divided by the number of processors.
  !! @param[in] matrix_dim the dimension of the actual matrix.
  !! @return a new dimension which includes padding.
  !! @todo write a more optimal algorithm using either plain brute force,
  !! or a prime factor based algorithm.
  PURE FUNCTION CalculateScaledDimension(this, matrix_dim) &
       & RESULT(scaled_dim)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: matrix_dim
    INTEGER :: scaled_dim
    !! Local Data
    INTEGER :: size_ratio
    INTEGER :: lcm

    lcm = this%process_grid%block_multiplier* &
         & this%process_grid%num_process_slices* &
         & this%process_grid%num_process_columns* &
         & this%process_grid%num_process_rows

    size_ratio = matrix_dim/lcm
    IF (size_ratio * lcm .EQ. matrix_dim) THEN
       scaled_dim = matrix_dim
    ELSE
       scaled_dim = (size_ratio + 1)*(lcm)
    END IF
  END FUNCTION CalculateScaledDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[in] matrix_to_split the matrix to split up.
  SUBROUTINE SplitMatrixToLocalBlocks_psr(this, matrix_to_split)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_lsr), INTENT(IN) :: matrix_to_split

#define LOCALMATRIX this%local_data_r
#include "includes/SplitMatrixToLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE SplitMatrixToLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[in] matrix_to_split the matrix to split up.
  SUBROUTINE SplitMatrixToLocalBlocks_psc(this, matrix_to_split)
    !! Parameters
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    TYPE(Matrix_lsc), INTENT(IN) :: matrix_to_split

#define LOCALMATRIX this%local_data_c
#include "includes/SplitMatrixToLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE SplitMatrixToLocalBlocks_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge together the local matrix blocks into one big matrix.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[inout] merged_matrix the merged matrix.
  PURE SUBROUTINE MergeMatrixLocalBlocks_psr(this, merged_matrix)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(Matrix_lsr), INTENT(INOUT) :: merged_matrix

#define LOCALMATRIX this%local_data_r
#include "includes/MergeMatrixLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE MergeMatrixLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge together the local matrix blocks into one big matrix.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[inout] merged_matrix the merged matrix.
  PURE SUBROUTINE MergeMatrixLocalBlocks_psc(this, merged_matrix)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN) :: this
    TYPE(Matrix_lsc), INTENT(INOUT) :: merged_matrix

#define LOCALMATRIX this%local_data_c
#include "includes/MergeMatrixLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE MergeMatrixLocalBlocks_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixPSModule
