!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Operations.
!! Unlike in earlier versions, this one will be blocked. Each Process
!! will hold many distributed sparse matrix blocks.
MODULE DistributedBlockedSparseMatrixModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DistributedMatrixMemoryPoolModule, ONLY : &
       & DistributedMatrixMemoryPool_t, ConstructDistributedMatrixMemoryPool, &
       & CheckDistributedMemoryPoolValidity
  USE GemmTasksModule
  USE MatrixGatherModule, ONLY : GatherHelper_t, GatherSizes, &
       & GatherAndComposeData, GatherAndComposeCleanup, GatherAndSumData, &
       & GatherAndSumCleanup, TestSizeRequest, TestOuterRequest, &
       & TestInnerRequest, TestDataRequest
  USE MatrixMarketModule
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation
  USE ProcessGridModule, ONLY : grid_error, IsRoot, RootID, total_processors,&
       & num_process_slices, num_process_rows, num_process_columns, &
       & slice_size, my_slice, my_row, my_column, global_rank,&
       & within_slice_rank, between_slice_rank, row_rank, column_rank, &
       & global_comm, within_slice_comm, row_comm, column_comm, &
       & between_slice_comm, &
       & blocked_column_comm, blocked_row_comm, blocked_between_slice_comm, &
       & number_of_blocks_rows, number_of_blocks_columns, block_multiplier
  USE SparseMatrixModule, ONLY : SparseMatrix_t, &
       & ConstructFromTripletList, DestructSparseMatrix, &
       & ConstructEmptySparseMatrix, CopySparseMatrix, &
       & TransposeSparseMatrix, SplitSparseMatrixColumns, &
       & ComposeSparseMatrixColumns, PrintSparseMatrix, &
       & MatrixToTripletList, DotSparseMatrix, PairwiseMultiplySparseMatrix, &
       & SparseMatrixNorm, ScaleSparseMatrix, IncrementSparseMatrix, Gemm, &
       & SparseMatrixGrandSum
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletModule, ONLY : Triplet_t, GetMPITripletType
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
       & DestructTripletList, SortTripletList, AppendToTripletList, &
       & SymmetrizeTripletList, GetTripletAt
  USE iso_c_binding
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for a distributed blocked CSR matrix.
  TYPE, PUBLIC :: DistributedSparseMatrix
     !> Number of matrix rows/columns for full matrix, scaled for process grid.
     INTEGER :: logical_matrix_dimension
     !> Number of matrix rows/columns for the full matrix, unscaled.
     INTEGER :: actual_matrix_dimension
     !! Local Storage
     !> A 2D array of local CSR matrices.
     TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: local_data
     INTEGER :: start_column !< first column stored locally.
     INTEGER :: end_column !< last column stored locally  is less than this.
     INTEGER :: start_row !< first row stored locally.
     INTEGER :: end_row !< last row stored locally is less than this.
     INTEGER :: local_columns !< number of local columns.
     INTEGER :: local_rows !< number of local rows.
  END TYPE DistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructors/Destructors
  PUBLIC :: ConstructEmpty
  PUBLIC :: DestructDistributedSparseMatrix
  PUBLIC :: CopyDistributedSparseMatrix
  !! File I/O
  PUBLIC :: ConstructFromMatrixMarket
  PUBLIC :: ConstructFromBinary
  !public :: ConstructFromDenseBinary
  PUBLIC :: WriteToMatrixMarket
  PUBLIC :: WriteToBinary
  !! Fill In Special Matrices
  PUBLIC :: FillFromTripletList
  PUBLIC :: FillDistributedIdentity
  PUBLIC :: FillDistributedPermutation
  !! Basic Accessors
  PUBLIC :: GetActualDimension
  PUBLIC :: GetLogicalDimension
  PUBLIC :: GetTripletList
  !! Basic Linear Algebra
  PUBLIC :: DotDistributedSparseMatrix
  PUBLIC :: IncrementDistributedSparseMatrix
  PUBLIC :: DistributedPairwiseMultiply
  PUBLIC :: DistributedGemm
  PUBLIC :: ScaleDistributedSparseMatrix
  PUBLIC :: DistributedSparseNorm
  PUBLIC :: Trace
  PUBLIC :: DistributedGrandSum
  PUBLIC :: EigenCircle
  PUBLIC :: ComputeSigma
  !! Utilities
  PUBLIC :: PrintDistributedSparseMatrix
  PUBLIC :: FilterSparseMatrix
  PUBLIC :: GetSize
  PUBLIC :: GetLoadBalance
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix.
  !! @param[out] this the matrix to be constructed.
  !! @param[in] matrix_dim_ the dimension of the full matrix.
  PURE SUBROUTINE ConstructEmpty(this, matrix_dim_)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    INTEGER, INTENT(in)           :: matrix_dim_

    CALL DestructDistributedSparseMatrix(this)

    !! Matrix Dimensions
    this%actual_matrix_dimension = matrix_dim_
    this%logical_matrix_dimension = CalculateScaledDimension(matrix_dim_, &
         & block_multiplier)
    !! Full Local Data Size Description
    this%local_rows = this%logical_matrix_dimension/num_process_rows
    this%local_columns = this%logical_matrix_dimension/num_process_columns
    !! Which Block Does This Process Hold?
    this%start_row = this%local_rows * my_row + 1
    this%end_row   = this%start_row + this%local_rows
    this%start_column = this%local_columns * my_column + 1
    this%end_column   = this%start_column + this%local_columns

    ALLOCATE(this%local_data(number_of_blocks_columns,number_of_blocks_rows))
  END SUBROUTINE ConstructEmpty
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix
  !! @param[in,out] this the matrix to destruct
  PURE SUBROUTINE DestructDistributedSparseMatrix(this)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    !! Local Data
    INTEGER :: row_counter, column_counter

    IF (ALLOCATED(this%local_data)) THEN
       DO row_counter = 1, number_of_blocks_rows
          DO column_counter = 1, number_of_blocks_columns
             CALL DestructSparseMatrix( &
                  & this%local_data(column_counter,row_counter))
          END DO
       END DO
       DEALLOCATE(this%local_data)
    END IF
  END SUBROUTINE DestructDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopyDistributedSparseMatrix(matA, matB)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)    :: matA
    TYPE(DistributedSparseMatrix), INTENT(inout) :: matB

    CALL DestructDistributedSparseMatrix(matB)
    matB = matA
  END SUBROUTINE CopyDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructFromMatrixMarket(this, file_name)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), INTENT(in) :: file_name
    INTEGER, PARAMETER :: MAX_LINE_LENGTH = 100
    !! File Handles
    INTEGER :: local_file_handler
    INTEGER :: mpi_file_handler
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Reading The File
    TYPE(TripletList_t) :: triplet_list
    TYPE(Triplet_t) :: temp_triplet
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

    !! Setup Involves Just The Root Opening And Reading Parameter Data
    CALL StartTimer("MPI Read Text")
    bytes_per_character = sizeof(temp_char)
    IF (IsRoot()) THEN
       header_length = 0
       local_file_handler = 16
       OPEN(local_file_handler, file=file_name, status="old")
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
    END IF

    !! Broadcast Parameters
    CALL MPI_Bcast(matrix_rows,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(matrix_columns,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(total_values,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(header_length,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(sparsity_type,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(data_type,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(pattern_type,1,MPI_INT,RootID,global_comm,grid_error)

    !! Build Local Storage
    CALL ConstructEmpty(this,matrix_rows)

    !! Global read
    CALL MPI_File_open(global_comm,file_name,MPI_MODE_RDONLY,&
         & MPI_INFO_NULL,mpi_file_handler,grid_error)

    !! Compute Offsets
    CALL MPI_File_get_size(mpi_file_handler,total_file_size,grid_error)
    local_data_size = (total_file_size - bytes_per_character*header_length)/&
         total_processors
    local_offset = bytes_per_character*header_length + &
         local_data_size*global_rank
    local_data_size_plus_buffer = local_data_size
    IF (.NOT. global_rank .EQ. total_processors-1) THEN
       local_data_size_plus_buffer = local_data_size + &
            MAX_LINE_LENGTH*bytes_per_character
    ELSE
       local_data_size_plus_buffer = (total_file_size - local_offset)
    END IF
    ALLOCATE(CHARACTER(LEN=local_data_size_plus_buffer) :: mpi_input_buffer)

    !! Do Actual Reading
    CALL MPI_File_read_at_all(mpi_file_handler,local_offset,mpi_input_buffer, &
         & INT(local_data_size_plus_buffer),MPI_CHARACTER,mpi_status,grid_error)

    !! Trim Off The Half Read Line At The Start
    IF (.NOT. global_rank .EQ. RootID) THEN
       full_buffer_counter = INDEX(mpi_input_buffer,new_line('A')) + 1
    ELSE
       full_buffer_counter = 1
    END IF

    !! Read By Line
    end_of_buffer = .FALSE.
    CALL ConstructTripletList(triplet_list)
    DO WHILE(.NOT. end_of_buffer)
       current_line_length = INDEX(mpi_input_buffer(full_buffer_counter:),&
            new_line('A'))

       IF (current_line_length .EQ. 0) THEN !! Hit The End Of The Buffer
          end_of_buffer = .TRUE.
       ELSE
          temp_substring = mpi_input_buffer( &
               & full_buffer_counter:full_buffer_counter+current_line_length-1)
          READ(temp_substring(:current_line_length-1),*) &
               & temp_triplet%index_row, temp_triplet%index_column, &
               & temp_triplet%point_value
          CALL AppendToTripletList(triplet_list, temp_triplet)

          IF (full_buffer_counter + current_line_length .GE. &
               & local_data_size+2) THEN
             IF (.NOT. global_rank .EQ. total_processors-1) THEN
                end_of_buffer = .TRUE.
             END IF
          END IF
          full_buffer_counter = full_buffer_counter + current_line_length
       END IF
    END DO

    !! Cleanup
    CALL MPI_File_close(mpi_file_handler,grid_error)
    CALL StopTimer("MPI Read Text")
    CALL MPI_Barrier(global_comm,grid_error)

    !! Redistribute The Matrix
    CALL SymmetrizeTripletList(triplet_list, pattern_type)
    CALL FillFromTripletList(this,triplet_list)

    CALL DestructTripletList(triplet_list)
    DEALLOCATE(mpi_input_buffer)
  END SUBROUTINE ConstructFromMatrixMarket
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructFromBinary(this, file_name)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), INTENT(in) :: file_name
    !! File Handles
    INTEGER :: mpi_file_handler
    !! Reading The File
    TYPE(TripletList_t) :: triplet_list
    INTEGER :: matrix_rows, matrix_columns, total_values
    INTEGER, DIMENSION(3) :: matrix_information
    INTEGER :: local_triplets
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER :: bytes_per_int, bytes_per_double
    INTEGER :: triplet_mpi_type
    !! Temporary variables
    INTEGER :: mpi_status(MPI_STATUS_SIZE)

    CALL StartTimer("MPI Read Binary")
    CALL MPI_Type_extent(MPI_INT,bytes_per_int,grid_error)
    CALL MPI_Type_extent(MPINTREAL,bytes_per_double,grid_error)

    CALL MPI_File_open(global_comm,file_name,MPI_MODE_RDONLY, &
         & MPI_INFO_NULL,mpi_file_handler,grid_error)

    !! Get The Matrix Parameters
    IF (IsRoot()) THEN
       local_offset = 0
       CALL MPI_File_read_at(mpi_file_handler, local_offset, &
            & matrix_information, 3, MPI_INT, mpi_status, grid_error)
       matrix_rows = matrix_information(1)
       matrix_columns = matrix_information(2)
       total_values = matrix_information(3)
    END IF

    !! Broadcast Parameters
    CALL MPI_Bcast(matrix_rows,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(matrix_columns,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(total_values,1,MPI_INT,RootID,global_comm,grid_error)

    !! Build Local Storage
    CALL ConstructEmpty(this, matrix_rows)

    !! Compute Offset
    local_triplets = total_values/total_processors
    local_offset = local_triplets * (global_rank)
    header_size = 3 * bytes_per_int
    IF (global_rank .EQ. total_processors - 1) THEN
       local_triplets = INT(total_values) - INT(local_offset)
    END IF
    local_offset = local_offset*(bytes_per_int*2+bytes_per_double) + header_size
    CALL ConstructTripletList(triplet_list,local_triplets)

    !! Do The Actual Reading
    triplet_mpi_type = GetMPITripletType()
    CALL MPI_File_set_view(mpi_file_handler,local_offset,triplet_mpi_type,&
         & triplet_mpi_type,"native",MPI_INFO_NULL,grid_error)
    CALL MPI_File_read_all(mpi_file_handler,triplet_list%data,local_triplets,&
         & triplet_mpi_type,mpi_status,grid_error)
    CALL MPI_File_close(mpi_file_handler,grid_error)
    CALL StopTimer("MPI Read Binary")

    CALL FillFromTripletList(this,triplet_list)

    !! Cleanup
    CALL DestructTripletList(triplet_list)
  END SUBROUTINE ConstructFromBinary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a file.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteToBinary(this,file_name)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), INTENT(in) :: file_name
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
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
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Determine Write Location
    bytes_per_int = sizeof(temp_int)
    bytes_per_double = sizeof(temp_double)
    local_data_size = SIZE(merged_local_data%values)*(bytes_per_int*2 + &
         & bytes_per_double*1)
    header_size = bytes_per_int*3
    ALLOCATE(local_values_buffer(slice_size))
    CALL MPI_Allgather(SIZE(merged_local_data%values),1,MPI_INT,&
         & local_values_buffer,1,MPI_INT,within_slice_comm,grid_error)
    write_offset = 0
    write_offset = write_offset + header_size
    DO counter = 1,within_slice_rank
       write_offset = write_offset + &
            & local_values_buffer(counter)*(bytes_per_int*2+bytes_per_double*1)
    END DO

    !! Write The File
    IF (between_slice_rank .EQ. 0) THEN
       !! Create Special MPI Type
       triplet_mpi_type = GetMPITripletType()

       CALL MatrixToTripletList(merged_local_data, triplet_list)
       !! Absolute Positions
       DO counter = 1, triplet_list%CurrentSize
          triplet_list%data(counter)%index_row = &
               & triplet_list%data(counter)%index_row + this%start_row - 1
          triplet_list%data(counter)%index_column = &
               & triplet_list%data(counter)%index_column + this%start_column - 1
       END DO
       CALL MPI_File_open(within_slice_comm,file_name,&
            & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL, &
            & mpi_file_handler, grid_error)
       !! Write Header
       IF (within_slice_rank .EQ. 0) THEN
          header_buffer(1) = this%actual_matrix_dimension
          header_buffer(2) = this%actual_matrix_dimension
          header_buffer(3) = SUM(local_values_buffer)
          CALL MPI_File_write_at(mpi_file_handler,zero_offset,header_buffer,3,&
               & MPI_INT,mpi_status,grid_error)
       END IF
       !! Write The Rest
       CALL MPI_File_set_view(mpi_file_handler,write_offset,triplet_mpi_type,&
            & triplet_mpi_type,"native",MPI_INFO_NULL,grid_error)
       CALL MPI_File_write(mpi_file_handler,triplet_list%data, &
            & triplet_list%CurrentSize, triplet_mpi_type,MPI_STATUS_IGNORE, &
            & grid_error)

       !! Cleanup
       CALL MPI_File_close(mpi_file_handler,grid_error)
       CALL DestructTripletList(triplet_list)
    END IF
    DEALLOCATE(local_values_buffer)
    CALL MPI_Barrier(global_comm,grid_error)
    CALL DestructSparseMatrix(merged_local_data)
  END SUBROUTINE WriteToBinary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a distributed sparse matrix to a matrix market file.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteToMatrixMarket(this,file_name)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), INTENT(in) :: file_name
    INTEGER, PARAMETER :: MAX_LINE_LENGTH = 1024
    !! Local MPI Variables
    INTEGER :: mpi_file_handler
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
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
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    bytes_per_character = sizeof(temp_char)

    !! Create the matrix size line
    NEW_LINE_LENGTH = LEN(new_line('A'))
    WRITE(temp_string1,'(A)') "%%MatrixMarket matrix coordinate real general" &
         & //new_line('A')//"%"//new_line('A')
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string1)) :: header_line1)
    header_line1 = TRIM(temp_string1)

    WRITE(temp_string2,*) this%actual_matrix_dimension, &
         & this%actual_matrix_dimension, GetSize(this)
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string2)+NEW_LINE_LENGTH) :: header_line2)
    WRITE(header_line2,*) TRIM(temp_string2)
    header_line2 = header_line2//new_line('A')

    header_size = LEN(header_line1) + LEN(header_line2)

    !! Local Data
    CALL MatrixToTripletList(merged_local_data, triplet_list)

    !! Absolute Positions
    DO counter = 1, triplet_list%CurrentSize
       triplet_list%data(counter)%index_row = &
            & triplet_list%data(counter)%index_row + this%start_row - 1
       triplet_list%data(counter)%index_column = &
            & triplet_list%data(counter)%index_column + this%start_column - 1
    END DO

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
       WRITE(write_buffer(offset_counter:offset_counter+temp_length),*) temp_string2(1:temp_length)
       offset_counter = offset_counter + temp_length
    END DO

    !! Figure out the offset sizes
    ALLOCATE(local_values_buffer(slice_size))
    CALL MPI_Allgather(triplet_list_string_length,1,MPI_INT,&
         & local_values_buffer,1,MPI_INT,within_slice_comm,grid_error)
    write_offset = 0
    write_offset = write_offset + header_size
    DO counter = 1,within_slice_rank
       write_offset = write_offset + &
            & local_values_buffer(counter)
    END DO

    !! Global Write
    IF (between_slice_rank .EQ. 0) THEN
       CALL MPI_File_open(within_slice_comm,file_name, &
            & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL, &
            & mpi_file_handler,grid_error)
       CALL MPI_File_set_size(mpi_file_handler,zero_size,grid_error)
       !! Write Header
       IF (within_slice_rank .EQ. 0) THEN
          header_offset = 0
          CALL MPI_File_write_at(mpi_file_handler,header_offset,header_line1, &
               & LEN(header_line1), MPI_CHARACTER,mpi_status,grid_error)
          header_offset = header_offset + LEN(header_line1)
          CALL MPI_File_write_at(mpi_file_handler,header_offset,header_line2, &
               & LEN(header_line2), MPI_CHARACTER,mpi_status,grid_error)
       END IF
       !! Write Local Data
       CALL MPI_File_set_view(mpi_file_handler,write_offset,MPI_CHARACTER,&
            & MPI_CHARACTER,"native",MPI_INFO_NULL,grid_error)
       CALL MPI_File_write(mpi_file_handler,write_buffer, &
            & triplet_list_string_length,&
            & MPI_CHARACTER,MPI_STATUS_IGNORE,grid_error)

       !! Cleanup
       CALL MPI_File_close(mpi_file_handler,grid_error)
    END IF
    CALL MPI_Barrier(global_comm,grid_error)
  END SUBROUTINE WriteToMatrixMarket
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !! should pass in triplet lists with global coordinates. It doesn't matter
  !! where each triplet is stored, as long as global coordinates are given.
  !! @param[inout] this the matrix to fill.
  !! @param[in] triplet_list the triplet list of values.
  SUBROUTINE FillFromTripletList(this,triplet_list)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    TYPE(TripletList_t) :: triplet_list
    !! Local Data
    TYPE(Permutation_t) :: basic_permutation
    TYPE(TripletList_t) :: sorted_triplet_list
    TYPE(SparseMatrix_t) :: local_matrix
    TYPE(SparseMatrix_t) :: gathered_matrix
    TYPE(GatherHelper_t) :: gather_helper
    REAL(ntreal), PARAMETER :: threshold = 0.0

    CALL StartTimer("FillFromTriplet")
    !! First we redistribute the triplet list to get all the local data
    !! on the correct process.
    CALL ConstructDefaultPermutation(basic_permutation, &
         & this%logical_matrix_dimension)
    CALL RedistributeTripletList(this,basic_permutation%index_lookup, &
         & basic_permutation%reverse_index_lookup, triplet_list, &
         & sorted_triplet_list)

    !! Now we can just construct a local matrix.
    CALL ConstructFromTripletList(local_matrix,sorted_triplet_list, &
         & this%local_rows, this%local_columns)

    !! And reduce over the Z dimension. This can be accomplished by
    !! summing up.
    CALL GatherSizes(local_matrix,between_slice_comm,gather_helper)
    DO WHILE(.NOT. TestSizeRequest(gather_helper))
    END DO
    CALL GatherAndSumData(local_matrix,between_slice_comm,gather_helper)
    DO WHILE(.NOT. TestOuterRequest(gather_helper))
    END DO
    DO WHILE(.NOT. TestInnerRequest(gather_helper))
    END DO
    DO WHILE(.NOT. TestDataRequest(gather_helper))
    END DO
    CALL GatherAndSumCleanUp(local_matrix, gathered_matrix, threshold, &
         & gather_helper)

    CALL SplitToLocalBlocks(this, gathered_matrix)
    CALL StopTimer("FillFromTriplet")

    CALL DestructSparseMatrix(local_matrix)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE FillFromTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] this the matrix being filled.
  PURE SUBROUTINE FillDistributedIdentity(this)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: unsorted_triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    INTEGER :: i, j
    INTEGER :: total_values
    TYPE(SparseMatrix_t) :: local_matrix

    !! There can't be more than one entry per row
    CALL ConstructTripletList(triplet_list,this%local_rows)

    total_values = 0
    !! Find local identity values
    row_iter: DO j = 1, this%local_rows
       column_iter: DO i = 1, this%local_columns
          IF (j + this%start_row - 1 .EQ. i + this%start_column - 1 .AND. &
               & j+this%start_row-1 .LE. this%actual_matrix_dimension) THEN
             total_values = total_values + 1
             triplet_list%data(total_values)%index_column = i
             triplet_list%data(total_values)%index_row = j
             triplet_list%data(total_values)%point_value = 1.0
          END IF
       END DO column_iter
    END DO row_iter

    !! Finish constructing
    CALL ConstructTripletList(unsorted_triplet_list,total_values)
    unsorted_triplet_list%data = triplet_list%data(:total_values)
    CALL SortTripletList(unsorted_triplet_list,this%local_columns,&
         & sorted_triplet_list)
    CALL ConstructFromTripletList(local_matrix,sorted_triplet_list, &
         & this%local_rows,this%local_columns)

    CALL SplitToLocalBlocks(this, local_matrix)

    CALL DestructSparseMatrix(local_matrix)
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(unsorted_triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE FillDistributedIdentity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !! If you don't specify permuterows, will default to permuting rows.
  !! @param[inout] this the matrix being filled.
  !! @param[in] permutation_vector describes for each row/column, where it goes.
  !! @param[in] permuterows if true permute rows, false permute columns.
  PURE SUBROUTINE FillDistributedPermutation(this, permutation_vector, &
       & permuterows)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    INTEGER, DIMENSION(:), INTENT(in) :: permutation_vector
    LOGICAL, OPTIONAL, INTENT(in) :: permuterows
    !! Local Data
    LOGICAL :: rows
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: unsorted_triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    INTEGER :: total_values
    INTEGER :: counter
    INTEGER :: local_row, local_column
    TYPE(SparseMatrix_t) :: local_matrix

    !! Figure out what type of permutation
    IF (PRESENT(permuterows) .AND. permuterows .EQV. .FALSE.) THEN
       rows = .FALSE.
    ELSE
       rows = .TRUE.
    END IF

    !! Build Local Triplet List
    !! There can't be more than one entry per row
    CALL ConstructTripletList(triplet_list,this%local_rows)
    total_values = 0
    IF (rows) THEN
       DO counter=this%start_row,this%end_row-1
          IF (permutation_vector(counter) .GE. this%start_column .AND. &
               & permutation_vector(counter) .LT. this%end_column) THEN
             total_values = total_values + 1
             local_column = permutation_vector(counter) - this%start_column + 1
             local_row = counter - this%start_row + 1
             triplet_list%data(total_values)%index_column = local_column
             triplet_list%data(total_values)%index_row = local_row
             triplet_list%data(total_values)%point_value = 1.0
          END IF
       END DO
    ELSE
       DO counter=this%start_column,this%end_column-1
          IF (permutation_vector(counter) .GE. this%start_row .AND. &
               & permutation_vector(counter) .LT. this%end_row) THEN
             total_values = total_values + 1
             local_column = counter - this%start_column + 1
             local_row = permutation_vector(counter) - this%start_row + 1
             triplet_list%data(total_values)%index_column = local_column
             triplet_list%data(total_values)%index_row = local_row
             triplet_list%data(total_values)%point_value = 1.0
          END IF
       END DO
    END IF

    !! Finish constructing
    CALL ConstructTripletList(unsorted_triplet_list,total_values)
    unsorted_triplet_list%data = triplet_list%data(:total_values)
    CALL SortTripletList(unsorted_triplet_list,this%local_columns,&
         & sorted_triplet_list)
    CALL ConstructFromTripletList(local_matrix,sorted_triplet_list, &
         & this%local_rows,this%local_columns)

    CALL SplitToLocalBlocks(this, local_matrix)

    CALL DestructSparseMatrix(local_matrix)
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(unsorted_triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE FillDistributedPermutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  PURE SUBROUTINE GetTripletList(this,triplet_list)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    TYPE(TripletList_t), INTENT(inout) :: triplet_list
    !! Local Data
    TYPE(SparseMatrix_t) :: merged_local_data
    INTEGER :: counter

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    CALL MatrixToTripletList(merged_local_data, triplet_list)
    DO counter=1, triplet_list%CurrentSize
       triplet_list%data(counter)%index_column = &
            & triplet_list%data(counter)%index_column + this%start_column - 1
       triplet_list%data(counter)%index_row = &
            & triplet_list%data(counter)%index_row + this%start_row - 1
    END DO
  END SUBROUTINE GetTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the actual dimension of the matrix.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetActualDimension(this) RESULT(DIMENSION)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%actual_matrix_dimension
  END FUNCTION GetActualDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the logical dimension of the matrix.
  !! @param[in] this the matrix.
  !! @result dimension of the matrix;
  PURE FUNCTION GetLogicalDimension(this) RESULT(DIMENSION)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    INTEGER :: DIMENSION
    DIMENSION = this%logical_matrix_dimension
  END FUNCTION GetLogicalDimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !! This will utilize the sparse vector increment routine.
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @param[in] alpha_in multiplier. Default value is 1.0
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0
  SUBROUTINE IncrementDistributedSparseMatrix(matA, matB, alpha_in,threshold_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix), INTENT(inout)  :: matB
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    !! Local Data
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: row_counter, column_counter

    !! Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0
    ELSE
       threshold = threshold_in
    END IF

    !$omp parallel
    !$omp do collapse(2)
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL IncrementSparseMatrix(matA%local_data(column_counter,row_counter),&
               & matB%local_data(column_counter,row_counter), alpha, threshold)
       END DO
    END DO
    !$omp end do
    !$omp end parallel

  END SUBROUTINE IncrementDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !! @param[in] matA Matrix A.
  !! @param[in,out] matB Matrix B.
  !! @result product the dot product.
  FUNCTION DotDistributedSparseMatrix(matA, matB) RESULT(product)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matB
    REAL(NTREAL) :: product
    !! Local Data
    TYPE(DistributedSparseMatrix)  :: matC

    CALL DistributedPairwiseMultiply(matA,matB,matC)
    product = DistributedGrandSum(matC)
    CALL DestructDistributedSparseMatrix(matC)
  END FUNCTION DotDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix.
  !! @param[in] matA Matrix A.
  !! @result sum the sum of all elements.
  FUNCTION DistributedGrandSum(matA) RESULT(sum)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matA
    REAL(NTREAL) :: sum
    !! Local Data
    INTEGER :: row_counter, column_counter
    REAL(NTREAL) :: temp

    sum = 0
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          temp = SparseMatrixGrandSum(matA%local_data(column_counter,row_counter))
          sum = sum + temp
       END DO
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
         & MPI_SUM, within_slice_comm, grid_error)
  END FUNCTION DistributedGrandSum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !! @param[in] matA Matrix A.
  !! @param[in] matB Matrix B.
  !! @param[in,out] matC = MatA mult MatB.
  SUBROUTINE DistributedPairwiseMultiply(matA, matB, matC)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matA
    TYPE(DistributedSparseMatrix), INTENT(in)  :: matB
    TYPE(DistributedSparseMatrix), INTENT(inout)  :: matC
    !! Local Data
    INTEGER :: row_counter, column_counter

    CALL ConstructEmpty(matC,matA%actual_matrix_dimension)

    !$omp parallel
    !$omp do collapse(2)
    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL PairwiseMultiplySparseMatrix( &
               & matA%local_data(column_counter,row_counter), &
               & matB%local_data(column_counter,row_counter), &
               & matC%local_data(column_counter,row_counter))
       END DO
    END DO
    !$omp end do
    !$omp end parallel
  END SUBROUTINE DistributedPairwiseMultiply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*matB+ beta*matC
  !! @param[in] matA Matrix A
  !! @param[in] matB Matrix B
  !! @param[out] matC = alpha*matA*matB + beta*matC
  !! @param[in] alpha_in scales the multiplication
  !! @param[in] beta_in scales matrix we sum on to
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @param[inout] memory_pool_in a memory pool that can be used for the
  !! calculation.
  SUBROUTINE DistributedGemm(matA,matB,matC,alpha_in,beta_in,threshold_in, &
       & memory_pool_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)    :: matA
    TYPE(DistributedSparseMatrix), INTENT(in)    :: matB
    TYPE(DistributedSparseMatrix), INTENT(inout) :: matC
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    TYPE(DistributedMatrixMemoryPool_t), OPTIONAL, INTENT(inout) :: &
         & memory_pool_in
    !! Local Versions of Optional Parameter
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(DistributedSparseMatrix) :: matAB
    !! Temporary Matrices
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: LocalRowContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution
    TYPE(SparseMatrix_t), DIMENSION(:,:), ALLOCATABLE :: SliceContribution
    !! Communication Helpers
    TYPE(GatherHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(GatherHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(GatherHelper_t), DIMENSION(:,:), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: row_counter, inner_row_counter
    INTEGER :: column_counter, inner_column_counter
    INTEGER :: duplicate_start_column, duplicate_offset_column
    INTEGER :: duplicate_start_row, duplicate_offset_row
    REAL(NTREAL) :: working_threshold
    !! Scheduling the A work
    INTEGER, DIMENSION(:), ALLOCATABLE :: ATasks
    INTEGER :: ATasks_completed
    !! Scheduling the B work
    INTEGER, DIMENSION(:), ALLOCATABLE :: BTasks
    INTEGER :: BTasks_completed
    !! Scheduling the AB work
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ABTasks
    INTEGER :: ABTasks_completed

    CALL StartTimer("GEMM")

    !! Handle the optional parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(beta_in)) THEN
       beta = 0.0
    ELSE
       beta = beta_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF
    !! The threshhold needs to be smaller if we're doing a sliced version
    !! because you might flush a value that would be kept in the summed version.
    IF (num_process_slices .GT. 1) THEN
       working_threshold = threshold/(num_process_slices*1000)
    ELSE
       working_threshold = threshold
    END IF

    !! Construct The Temporary Matrices
    CALL ConstructEmpty(matAB,matA%actual_matrix_dimension)

    ALLOCATE(AdjacentABlocks(number_of_blocks_columns/num_process_slices, &
         & number_of_blocks_rows))
    ALLOCATE(LocalRowContribution(number_of_blocks_rows))
    ALLOCATE(GatheredRowContribution(number_of_blocks_rows))
    ALLOCATE(GatheredRowContributionT(number_of_blocks_rows))

    ALLOCATE(TransposedBBlocks(number_of_blocks_columns, &
         & number_of_blocks_rows/num_process_slices))
    ALLOCATE(LocalColumnContribution(number_of_blocks_columns))
    ALLOCATE(GatheredColumnContribution(number_of_blocks_columns))
    ALLOCATE(SliceContribution(number_of_blocks_columns, &
         & number_of_blocks_rows))

    !! Helpers
    ALLOCATE(row_helper(number_of_blocks_rows))
    ALLOCATE(column_helper(number_of_blocks_columns))
    ALLOCATE(slice_helper(number_of_blocks_columns, &
         & number_of_blocks_rows))

    !! Construct the task queues
    ALLOCATE(ATasks(number_of_blocks_rows))
    DO row_counter=1,number_of_blocks_rows
       ATasks(row_counter) = LocalGatherA
    END DO
    ALLOCATE(BTasks(number_of_blocks_columns))
    DO column_counter=1,number_of_blocks_columns
       BTasks(column_counter) = LocalGatherB
    END DO
    ALLOCATE(ABTasks(number_of_blocks_columns,number_of_blocks_rows))
    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          ABTasks(column_counter,row_counter) = AwaitingAB
       END DO
    END DO

    !! Setup A Tasks
    duplicate_start_column = my_slice+1
    duplicate_offset_column = num_process_slices

    !! Setup B Tasks
    duplicate_start_row = my_slice+1
    duplicate_offset_row = num_process_slices

    !! Setup AB Tasks
    IF (PRESENT(memory_pool_in)) THEN
       IF (.NOT. CheckDistributedMemoryPoolValidity(memory_pool_in)) THEN
          CALL ConstructDistributedMatrixMemoryPool(memory_pool_in)
       END IF
    END IF

    !! Run A Tasks
    ATasks_completed = 0
    BTasks_completed = 0
    ABTasks_completed = 0
    !$omp PARALLEL
    !$omp MASTER
    DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
         & BTasks_completed .LT. SIZE(BTasks) .OR. &
         & ABTasks_completed .LT. SIZE(ABTasks))
       DO row_counter=1,number_of_blocks_rows
          SELECT CASE (ATasks(row_counter))
          CASE(LocalGatherA)
             ATasks(row_counter) = TaskRunningA
             !$omp task default(shared), private(inner_column_counter), &
             !$omp& firstprivate(row_counter)
             !! First Align The Data We're Working With
             DO inner_column_counter=1, &
                  & number_of_blocks_columns/num_process_slices
                CALL CopySparseMatrix(matA%local_data(duplicate_start_column+ &
                     & duplicate_offset_column*(inner_column_counter-1), &
                     & row_counter), &
                     & AdjacentABlocks(inner_column_counter,row_counter))
             END DO
             !! Then Do A Local Gather
             CALL ComposeSparseMatrixColumns(AdjacentABlocks(:,row_counter), &
                  & LocalRowContribution(row_counter))
             ATasks(row_counter) = SendSizeA
             !$omp end task
          CASE(SendSizeA)
             !! Then Start A Global Gather
             CALL GatherSizes(LocalRowContribution(row_counter), &
                  & blocked_row_comm(row_counter), row_helper(row_counter))
             ATasks(row_counter) = ComposeA
          CASE(ComposeA)
             IF (TestSizeRequest(row_helper(row_counter))) THEN
                CALL GatherAndComposeData(LocalRowContribution(row_counter), &
                     & blocked_row_comm(row_counter), &
                     & GatheredRowContribution(row_counter), &
                     & row_helper(row_counter))
                ATasks(row_counter) = WaitOuterA
             END IF
          CASE(WaitOuterA)
             IF (TestOuterRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = WaitInnerA
             END IF
          CASE(WaitInnerA)
             IF (TestInnerRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = WaitDataA
             END IF
          CASE(WaitDataA)
             IF (TestDataRequest(row_helper(row_counter))) THEN
                ATasks(row_counter) = AdjustIndicesA
             END IF
          CASE(AdjustIndicesA)
             ATasks(row_counter) = TaskRunningA
             !$omp task default(shared), firstprivate(row_counter)
             CALL GatherAndComposeCleanup(LocalRowContribution(row_counter), &
                  & GatheredRowContribution(row_counter), &
                  & row_helper(row_counter))
             CALL TransposeSparseMatrix(GatheredRowContribution(row_counter), &
                  & GatheredRowContributionT(row_counter))
             ATasks(row_counter) = CleanupA
             !$omp end task
          CASE(CleanupA)
             ATasks(row_counter) = FinishedA
             ATasks_completed = ATasks_completed + 1
          END SELECT
       END DO
       !! B Tasks
       DO column_counter=1,number_of_blocks_columns
          SELECT CASE (BTasks(column_counter))
          CASE(LocalGatherB)
             BTasks(column_counter) = TaskRunningB
             !$omp task default(shared), private(inner_row_counter), &
             !$omp& firstprivate(column_counter)
             !! First Transpose The Data We're Working With
             DO inner_row_counter=1,number_of_blocks_rows/num_process_slices
                CALL TransposeSparseMatrix(matB%local_data(column_counter, &
                     & duplicate_start_row+ &
                     & duplicate_offset_row*(inner_row_counter-1)), &
                     & TransposedBBlocks(column_counter,inner_row_counter))
             END DO
             !! Then Do A Local Gather
             CALL ComposeSparseMatrixColumns(&
                  & TransposedBBlocks(column_counter,:),&
                  & LocalColumnContribution(column_counter))
             BTasks(column_counter) = SendSizeB
             !$omp end task
          CASE(SendSizeB)
             !! Then A Global Gather
             CALL GatherSizes(LocalColumnContribution(column_counter), &
                  & blocked_column_comm(column_counter), &
                  & column_helper(column_counter))
             BTasks(column_counter) = LocalComposeB
          CASE(LocalComposeB)
             IF (TestSizeRequest(column_helper(column_counter))) THEN
                CALL GatherAndComposeData( &
                     & LocalColumnContribution(column_counter),&
                     & blocked_column_comm(column_counter), &
                     & GatheredColumnContribution(column_counter), &
                     & column_helper(column_counter))
                BTasks(column_counter) = WaitOuterB
             END IF
          CASE(WaitOuterB)
             IF (TestOuterRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = WaitInnerB
             END IF
          CASE(WaitInnerB)
             IF (TestInnerRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = WaitDataB
             END IF
          CASE(WaitDataB)
             IF (TestDataRequest(column_helper(column_counter))) THEN
                BTasks(column_counter) = AdjustIndicesB
             END IF
          CASE(AdjustIndicesB)
             BTasks(column_counter) = TaskRunningB
             !$omp task default(shared), firstprivate(column_counter)
             CALL GatherAndComposeCleanup( &
                  & LocalColumnContribution(column_counter),&
                  & GatheredColumnContribution(column_counter), &
                  & column_helper(column_counter))
             BTasks(column_counter) = CleanupB
             !$omp end task
          CASE(CleanupB)
             BTasks(column_counter) = FinishedB
             BTasks_completed = BTasks_completed + 1
          END SELECT
       END DO
       !! AB Tasks
       DO row_counter=1,number_of_blocks_rows
          DO column_counter=1,number_of_blocks_columns
             SELECT CASE(ABTasks(column_counter,row_counter))
             CASE (AwaitingAB)
                IF (ATasks(row_counter) .EQ. FinishedA .AND. &
                     & BTasks(column_counter) .EQ. FinishedB) THEN
                   ABTasks(column_counter,row_counter) = GemmAB
                END IF
             CASE (GemmAB)
                ABTasks(column_counter,row_counter) = TaskRunningAB
                !$omp task default(shared), &
                !$omp& firstprivate(row_counter,column_counter)
                IF (PRESENT(memory_pool_in)) THEN
                   CALL Gemm(GatheredRowContributionT(row_counter), &
                        & GatheredColumnContribution(column_counter), &
                        & SliceContribution(column_counter,row_counter), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold, &
                        & blocked_memory_pool_in= &
                        & memory_pool_in%grid(column_counter,row_counter))
                ELSE
                   CALL Gemm(GatheredRowContributionT(row_counter), &
                        & GatheredColumnContribution(column_counter), &
                        & SliceContribution(column_counter,row_counter), &
                        & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                        & alpha_in=alpha, threshold_in=working_threshold)
                END IF
                ABTasks(column_counter,row_counter) = SendSizeAB
                !$omp end task
             CASE(SendSizeAB)
                CALL GatherSizes(SliceContribution(column_counter,row_counter), &
                     & blocked_between_slice_comm(column_counter,row_counter), &
                     & slice_helper(column_counter,row_counter))
                ABTasks(column_counter,row_counter) = GatherAndSumAB
             CASE (GatherAndSumAB)
                IF (TestSizeRequest(slice_helper(column_counter,row_counter))) THEN
                   CALL GatherAndSumData( &
                        & SliceContribution(column_counter,row_counter), &
                        & blocked_between_slice_comm(column_counter,row_counter), &
                        & slice_helper(column_counter,row_counter))
                   ABTasks(column_counter,row_counter) = WaitOuterAB
                END IF
             CASE (WaitOuterAB)
                IF (TestOuterRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = WaitInnerAB
                END IF
             CASE (WaitInnerAB)
                IF (TestInnerRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = WaitDataAB
                END IF
             CASE (WaitDataAB)
                IF (TestDataRequest(slice_helper(column_counter,row_counter))) THEN
                   ABTasks(column_counter,row_counter) = LocalSumAB
                END IF
             CASE(LocalSumAB)
                ABTasks(column_counter,row_counter) = TaskRunningAB
                !$omp task default(shared), &
                !$omp& firstprivate(row_counter,column_counter)
                CALL GatherAndSumCleanup( &
                     & SliceContribution(column_counter,row_counter), &
                     & matAB%local_data(column_counter,row_counter), &
                     & threshold, slice_helper(column_counter,row_counter))
                ABTasks(column_counter,row_counter) = CleanupAB
                !$omp end task
             CASE(CleanupAB)
                ABTasks(column_counter,row_counter) = FinishedAB
                ABTasks_completed = ABTasks_completed + 1
             END SELECT
          END DO
       END DO
    END DO
    !$omp end MASTER
    !$omp end PARALLEL

    !! Copy to output matrix.
    IF (beta .EQ. 0.0) THEN
       CALL CopyDistributedSparseMatrix(matAB,matC)
    ELSE
       CALL ScaleDistributedSparseMatrix(MatC,beta)
       CALL IncrementDistributedSparseMatrix(MatAB,MatC)
    END IF

    !! Cleanup
    CALL DestructDistributedSparseMatrix(matAB)
    DEALLOCATE(row_helper)
    DEALLOCATE(column_helper)
    DEALLOCATE(slice_helper)

    !! Deallocate Buffers From A
    DO row_counter=1,number_of_blocks_rows
       DO inner_column_counter=1,number_of_blocks_columns/num_process_slices
          CALL DestructSparseMatrix(AdjacentABlocks(inner_column_counter, &
               & row_counter))
       END DO
       CALL DestructSparseMatrix(LocalRowContribution(row_counter))
       CALL DestructSparseMatrix(GatheredRowContribution(row_counter))
    END DO
    DEALLOCATE(AdjacentABlocks)
    DEALLOCATE(LocalRowContribution)
    DEALLOCATE(GatheredRowContribution)
    !! Deallocate Buffers From B
    DO column_counter=1,number_of_blocks_columns
       DO inner_row_counter=1,number_of_blocks_rows/num_process_slices
          CALL DestructSparseMatrix(TransposedBBlocks(column_counter, &
               & inner_row_counter))
       END DO
       CALL DestructSparseMatrix(LocalColumnContribution(column_counter))
    END DO
    DEALLOCATE(TransposedBBlocks)
    DEALLOCATE(LocalColumnContribution)
    !! Deallocate Buffers From Multiplying The Block
    DO row_counter=1,number_of_blocks_rows
       CALL DestructSparseMatrix(GatheredRowContributionT(row_counter))
    END DO
    DO column_counter=1,number_of_blocks_columns
       CALL DestructSparseMatrix(GatheredColumnContribution(column_counter))
    END DO
    DEALLOCATE(GatheredRowContributionT)
    DEALLOCATE(GatheredColumnContribution)
    !! Deallocate Buffers From Sum
    DO row_counter=1,number_of_blocks_rows
       DO column_counter=1,number_of_blocks_columns
          CALL DestructSparseMatrix(SliceContribution(column_counter,row_counter))
       END DO
    END DO
    DEALLOCATE(SliceContribution)

    CALL StopTimer("GEMM")
  END SUBROUTINE DistributedGemm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  !! @param[inout] this Matrix to scale.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleDistributedSparseMatrix(this,constant)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    REAL(NTREAL), INTENT(in) :: constant
    INTEGER :: row_counter, column_counter

    DO row_counter = 1, number_of_blocks_rows
       DO column_counter = 1, number_of_blocks_columns
          CALL ScaleSparseMatrix(this%local_data(column_counter,row_counter), &
               & constant)
       END DO
    END DO
  END SUBROUTINE ScaleDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  !! @param[in] this the matrix to compute the norm of.
  !! @return the norm value of the full distributed sparse matrix.
  FUNCTION DistributedSparseNorm(this) RESULT(norm_value)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Sum Along Columns
    CALL SparseMatrixNorm(merged_local_data,local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
         & MPINTREAL, MPI_SUM,column_comm,grid_error)

    !! Find Max Value Amonst Columns
    norm_value = MAXVAL(local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
         & row_comm, grid_error)

    CALL DestructSparseMatrix(merged_local_data)
    DEALLOCATE(local_norm)
  END FUNCTION DistributedSparseNorm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  !! @param[in] this the matrix to compute the norm of.
  !! @return the trace value of the full distributed sparse matrix.
  FUNCTION Trace(this) RESULT(trace_value)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    REAL(NTREAL) :: trace_value
    !! Local data
    TYPE(TripletList_t) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Compute The Local Contribution
    trace_value = 0
    CALL MatrixToTripletList(merged_local_data,triplet_list)
    DO counter = 1, triplet_list%CurrentSize
       IF (this%start_row + triplet_list%data(counter)%index_row .EQ. &
            & this%start_column + triplet_list%data(counter)%index_column) THEN
          trace_value = trace_value + triplet_list%data(counter)%point_value
       END IF
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
         & MPI_SUM, within_slice_comm,grid_error)

    CALL DestructSparseMatrix(merged_local_data)
  END FUNCTION Trace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !! Uses Gershgorin's theorem.
  !! @param[in] this the matrix to compute the min/max of.
  !! @param[out] min_value a lower bound on the eigenspectrum.
  !! @param[out] max_value an uppder bound on the eigenspectrum.
  SUBROUTINE EigenCircle(this,min_value,max_value)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    REAL(NTREAL), INTENT(out) :: min_value, max_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    TYPE(TripletList_t) :: triplet_list
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Allocate Space For Result
    ALLOCATE(per_column_min(merged_local_data%columns))
    ALLOCATE(per_column_max(merged_local_data%columns))

    !! Compute The Local Contribution
    per_column_min = 0
    per_column_max = 0
    CALL MatrixToTripletList(merged_local_data,triplet_list)
    DO counter = 1, triplet_list%CurrentSize
       IF (this%start_row + triplet_list%data(counter)%index_row .EQ. &
            & this%start_column + triplet_list%data(counter)%index_column) THEN
          per_column_min(triplet_list%data(counter)%index_column) = &
               & per_column_min(triplet_list%data(counter)%index_column) + &
               & triplet_list%data(counter)%point_value
          per_column_max(triplet_list%data(counter)%index_column) = &
               & per_column_max(triplet_list%data(counter)%index_column) + &
               & triplet_list%data(counter)%point_value
       ELSE
          per_column_min(triplet_list%data(counter)%index_column) = &
               & per_column_min(triplet_list%data(counter)%index_column) - &
               & ABS(triplet_list%data(counter)%point_value)
          per_column_max(triplet_list%data(counter)%index_column) = &
               & per_column_max(triplet_list%data(counter)%index_column) + &
               & ABS(triplet_list%data(counter)%point_value)
       END IF
    END DO

    !! Sum Along Columns
    CALL MPI_Allreduce(MPI_IN_PLACE,per_column_min,SIZE(per_column_min), &
         & MPINTREAL,MPI_SUM,column_comm,grid_error)
    CALL MPI_Allreduce(MPI_IN_PLACE,per_column_max,SIZE(per_column_max), &
         & MPINTREAL,MPI_SUM,column_comm,grid_error)

    min_value = MINVAL(per_column_min)
    max_value = MAXVAL(per_column_max)

    CALL MPI_Allreduce(MPI_IN_PLACE,min_value,1,MPINTREAL,MPI_MIN, &
         & row_comm, grid_error)
    CALL MPI_Allreduce(MPI_IN_PLACE,max_value,1,MPINTREAL,MPI_MAX, &
         & row_comm, grid_error)

    DEALLOCATE(per_column_min)
    DEALLOCATE(per_column_max)
    CALL DestructSparseMatrix(merged_local_data)
  END SUBROUTINE EigenCircle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !! @todo describe this better.
  !! @param[in] this the matrix to compute the sigma value of.
  !! @param[out] sigma_value sigma.
  SUBROUTINE ComputeSigma(this,sigma_value)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    REAL(NTREAL), INTENT(out) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: &
         & column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    ALLOCATE(column_sigma_contribution(merged_local_data%columns))
    column_sigma_contribution = 0
    DO outer_counter = 1, merged_local_data%columns
       DO inner_counter = merged_local_data%outer_index(outer_counter), &
            & merged_local_data%outer_index(outer_counter+1)-1
          column_sigma_contribution(outer_counter) = &
               & column_sigma_contribution(outer_counter) + &
               & ABS(merged_local_data%values(inner_counter+1))
       END DO
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE,column_sigma_contribution,&
         & merged_local_data%columns,MPINTREAL,MPI_SUM, &
         & column_comm, grid_error)
    CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
         & MPINTREAL,MPI_MAX,row_comm,grid_error)
    sigma_value = 1.0d+0/(sigma_value**2)

    DEALLOCATE(column_sigma_contribution)
    CALL DestructSparseMatrix(merged_local_data)
  END SUBROUTINE ComputeSigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print ouf a distributed sparse matrix.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintDistributedSparseMatrix(this, file_name_in)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: file_name_in
    !! Helpers For Communication
    TYPE(GatherHelper_t) :: row_helper
    TYPE(GatherHelper_t) :: column_helper
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    !! Temporary Variables
    TYPE(SparseMatrix_t) :: merged_local_data
    TYPE(SparseMatrix_t) :: merged_local_dataT
    TYPE(SparseMatrix_t) :: merged_columns
    TYPE(SparseMatrix_t) :: merged_columnsT
    TYPE(SparseMatrix_t) :: full_gathered

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    !! Merge Columns
    CALL TransposeSparseMatrix(merged_local_data,merged_local_dataT)
    CALL GatherSizes(merged_local_dataT, column_comm, column_helper)
    CALL MPI_Wait(column_helper%size_request,mpi_status,grid_error)
    CALL GatherAndComposeData(merged_local_dataT,column_comm,merged_columns, &
         & column_helper)
    CALL MPI_Wait(column_helper%outer_request,mpi_status,grid_error)
    CALL MPI_Wait(column_helper%inner_request,mpi_status,grid_error)
    CALL MPI_Wait(column_helper%data_request,mpi_status,grid_error)
    CALL GatherAndComposeCleanup(merged_local_dataT,merged_columns, &
         & column_helper)

    !! Merge Rows
    CALL TransposeSparseMatrix(merged_columns,merged_columnsT)
    !call GatherAllAndCompose(merged_columnsT, row_comm, full_gathered)
    CALL GatherSizes(merged_columnsT, row_comm, row_helper)
    CALL MPI_Wait(row_helper%size_request,mpi_status,grid_error)
    CALL GatherAndComposeData(merged_columnsT,row_comm,full_gathered,row_helper)
    CALL MPI_Wait(row_helper%outer_request,mpi_status,grid_error)
    CALL MPI_Wait(row_helper%inner_request,mpi_status,grid_error)
    CALL MPI_Wait(row_helper%data_request,mpi_status,grid_error)
    CALL GatherAndComposeCleanup(merged_columnsT,full_gathered,row_helper)

    !! Make these changes so that it prints the logical rows/columns
    full_gathered%rows = this%actual_matrix_dimension
    full_gathered%columns = this%actual_matrix_dimension

    IF (IsRoot()) THEN
       IF (PRESENT(file_name_in)) THEN
          CALL PrintSparseMatrix(full_gathered, file_name_in)
       ELSE
          CALL PrintSparseMatrix(full_gathered)
       END IF
    END IF

    CALL DestructSparseMatrix(merged_local_data)
    CALL DestructSparseMatrix(merged_local_dataT)
    CALL DestructSparseMatrix(merged_columns)
    CALL DestructSparseMatrix(merged_columnsT)
  END SUBROUTINE PrintDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !! All (absolute) values below the threshold are set to zero.
  !! @param[inout] this matrix to filter
  !! @param[in] threshold (absolute) values below this are filtered
  SUBROUTINE FilterSparseMatrix(this, threshold)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: new_list
    TYPE(Triplet_t) :: temporary
    INTEGER :: counter
    INTEGER :: size_temp

    CALL GetTripletList(this,triplet_list)
    CALL ConstructTripletList(new_list)
    DO counter=1,triplet_list%CurrentSize
      CALL GetTripletAt(triplet_list,counter,temporary)
      IF (ABS(temporary%point_value) .GT. threshold) THEN
        CALL AppendToTripletList(new_list,temporary)
      END IF
    END DO
    size_temp = this%actual_matrix_dimension
    CALL DestructDistributedSparseMatrix(this)
    CALL ConstructEmpty(this,size_temp)
    CALL FillFromTripletList(this,new_list)
  END SUBROUTINE FilterSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the total number of non zero entries in the distributed sparse matrix.
  !! @param[in] this the distributed sparse matrix to calculate the non-zero
  !! entries of.
  !! @return the number of non-zero entries in the matrix.
  FUNCTION GetSize(this) RESULT(total_size)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    INTEGER(c_long) :: total_size
    !! Local Data
    !integer :: local_size
    REAL(NTREAL) :: local_size
    REAL(NTREAL) :: temp_size
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    local_size = SIZE(merged_local_data%values)
    CALL MPI_Allreduce(local_size,temp_size,1,MPINTREAL,MPI_SUM,&
         & within_slice_comm, grid_error)

    total_size = INT(temp_size,kind=c_long)

    CALL DestructSparseMatrix(merged_local_data)
  END FUNCTION GetSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get a measure of how load balanced this matrix is. For each process, the
  !! number of non-zero entries is calculated. Then, this function returns
  !! the max and min of those values.
  !! @param[in] this The matrix to compute the measure on.
  !! @param[out] min_size the minimum entries contained on a single process.
  !! @param[out] max_size the maximum entries contained on a single process.
  SUBROUTINE GetLoadBalance(this, min_size, max_size)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    INTEGER, INTENT(out) :: min_size
    INTEGER, INTENT(out) :: max_size
    !! Local Data
    INTEGER :: local_size
    TYPE(SparseMatrix_t) :: merged_local_data

    !! Merge all the local data
    CALL MergeLocalBlocks(this, merged_local_data)

    local_size = SIZE(merged_local_data%values)
    CALL MPI_Allreduce(local_size,max_size,1,MPI_INT,MPI_MAX,&
         & within_slice_comm, grid_error)
    CALL MPI_Allreduce(local_size,min_size,1,MPI_INT,MPI_MIN,&
         & within_slice_comm, grid_error)

    CALL DestructSparseMatrix(merged_local_data)
  END SUBROUTINE GetLoadBalance
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
  SUBROUTINE RedistributeTripletList(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    INTEGER, DIMENSION(:), INTENT(in) :: index_lookup
    INTEGER, DIMENSION(:), INTENT(in) :: reverse_index_lookup
    TYPE(TripletList_t), INTENT(in) :: initial_triplet_list
    TYPE(TripletList_t), INTENT(out) :: sorted_triplet_list
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_lookup
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_lookup
    INTEGER, DIMENSION(:), ALLOCATABLE :: location_list_within_slice
    INTEGER, DIMENSION(:), ALLOCATABLE :: internal_offset
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_offset
    INTEGER, DIMENSION(:), ALLOCATABLE :: receive_offset
    TYPE(TripletList_t) :: gathered_triplet_list
    !! Send Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_index_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_index_column
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_values
    INTEGER, DIMENSION(:), ALLOCATABLE :: receive_buffer_index_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: receive_buffer_index_column
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: receive_buffer_values
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_count
    INTEGER, DIMENSION(:), ALLOCATABLE :: receive_count
    !! Temporary Values
    INTEGER :: row_size, column_size
    INTEGER :: temp_row, temp_column
    INTEGER :: temp_index
    INTEGER :: process_id
    TYPE(Triplet_t) :: temp_triplet
    INTEGER :: counter

    CALL StartTimer("Redistribute")
    !! First we need to figure out where our local elements go
    ALLOCATE(row_lookup(SIZE(index_lookup)))
    ALLOCATE(column_lookup(SIZE(index_lookup)))
    row_size = SIZE(index_lookup)/num_process_rows
    DO counter = LBOUND(index_lookup,1), UBOUND(index_lookup,1)
       row_lookup(index_lookup(counter)) = (counter-1)/(row_size)
    END DO
    column_size = SIZE(index_lookup)/num_process_columns
    DO counter = LBOUND(index_lookup,1), UBOUND(index_lookup,1)
       column_lookup(index_lookup(counter)) = (counter-1)/(column_size)
    END DO
    ALLOCATE(location_list_within_slice(initial_triplet_list%CurrentSize))
    DO counter = 1, initial_triplet_list%CurrentSize
       temp_row = row_lookup(initial_triplet_list%data(counter)%index_row)
       temp_column = &
            & column_lookup(initial_triplet_list%data(counter)%index_column)
       location_list_within_slice(counter) = &
            & temp_column+temp_row*num_process_columns
    END DO

    !! To Each Process, How Many Items To Send, How Many To Receive
    ALLOCATE(send_count(slice_size))
    ALLOCATE(receive_count(slice_size))
    send_count = 0
    DO counter = 1, SIZE(location_list_within_slice)
       temp_index = location_list_within_slice(counter)
       send_count(temp_index+1) = send_count(temp_index+1) + 1
    END DO
    ALLOCATE(send_offset(slice_size+1))
    send_offset(1) = 0
    DO counter = 1,slice_size
       send_offset(counter+1) = send_offset(counter) + send_count(counter)
    END DO
    CALL MPI_Alltoall(send_count,1,MPI_INT,receive_count,1,MPI_INT,&
         & within_slice_comm,grid_error)
    ALLOCATE(receive_offset(slice_size+1))
    receive_offset(1) = 0
    DO counter = 1, slice_size
       receive_offset(counter+1)=receive_offset(counter)+receive_count(counter)
    END DO
    ALLOCATE(receive_buffer_index_row(SUM(receive_count)))
    ALLOCATE(receive_buffer_index_column(SUM(receive_count)))
    ALLOCATE(receive_buffer_values(SUM(receive_count)))

    !! Prepare the send buffers (and implicitly sort)
    ALLOCATE(send_buffer_index_row(initial_triplet_list%CurrentSize))
    ALLOCATE(send_buffer_index_column(initial_triplet_list%CurrentSize))
    ALLOCATE(send_buffer_values(initial_triplet_list%CurrentSize))
    ALLOCATE(internal_offset(slice_size))
    internal_offset = 0
    DO counter = 1, initial_triplet_list%CurrentSize
       process_id = location_list_within_slice(counter)
       temp_index = send_offset(process_id+1) + internal_offset(process_id+1) + 1
       send_buffer_index_row(temp_index) = &
            & initial_triplet_list%data(counter)%index_row
       send_buffer_index_column(temp_index) = &
            & initial_triplet_list%data(counter)%index_column
       send_buffer_values(temp_index) = &
            & initial_triplet_list%data(counter)%point_value
       internal_offset(process_id+1) = internal_offset(process_id+1) + 1
    END DO

    !! Send
    CALL MPI_Alltoallv(send_buffer_index_row,send_count,send_offset,MPI_INT,&
         & receive_buffer_index_row,receive_count,receive_offset,MPI_INT,&
         & within_slice_comm,grid_error)
    CALL MPI_Alltoallv(send_buffer_index_column,send_count,send_offset,MPI_INT,&
         & receive_buffer_index_column,receive_count,receive_offset,MPI_INT,&
         & within_slice_comm,grid_error)
    CALL MPI_Alltoallv(send_buffer_values,send_count,send_offset, &
         & MPINTREAL,receive_buffer_values,receive_count, &
         & receive_offset,MPINTREAL,within_slice_comm,grid_error)

    !! Construct Matrix
    CALL ConstructTripletList(gathered_triplet_list,SIZE(receive_buffer_values))
    DO counter = 1, SIZE(receive_buffer_values)
       temp_triplet%index_row = &
            & reverse_index_lookup(receive_buffer_index_row(counter)) - &
            & this%start_row + 1
       temp_triplet%index_column = &
            & reverse_index_lookup(receive_buffer_index_column(counter)) - &
            & this%start_column + 1
       temp_triplet%point_value = receive_buffer_values(counter)
       gathered_triplet_list%data(counter) = temp_triplet
    END DO
    CALL SortTripletList(gathered_triplet_list,this%local_columns,&
         & sorted_triplet_list)

    !! Cleanup
    DEALLOCATE(row_lookup)
    DEALLOCATE(column_lookup)
    DEALLOCATE(location_list_within_slice)
    DEALLOCATE(internal_offset)
    DEALLOCATE(send_offset)
    DEALLOCATE(receive_offset)
    CALL DestructTripletList(gathered_triplet_list)
    DEALLOCATE(send_buffer_index_row)
    DEALLOCATE(send_buffer_index_column)
    DEALLOCATE(send_buffer_values)
    DEALLOCATE(receive_buffer_index_row)
    DEALLOCATE(receive_buffer_index_column)
    DEALLOCATE(receive_buffer_values)
    DEALLOCATE(send_count)
    DEALLOCATE(receive_count)

    CALL StopTimer("Redistribute")
  END SUBROUTINE RedistributeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate a matrix size that can be divided by the number of processors.
  !! @param[in] matrix_dim the dimension of the actual matrix.
  !! @return a new dimension which includes padding.
  !! @todo write a more optimal algorithm using either plain brute force,
  !! or a prime factor based algorithm.
  PURE FUNCTION CalculateScaledDimension(matrix_dim, block_multiplier) &
       & RESULT(scaled_dim)
    !! Parameters
    INTEGER, INTENT(in) :: matrix_dim
    INTEGER, INTENT(in) :: block_multiplier
    INTEGER :: scaled_dim
    !! Local Data
    INTEGER :: size_ratio
    INTEGER :: lcm

    lcm = block_multiplier*num_process_slices* &
         & num_process_columns*num_process_rows
    !lcm = ComputeLCM(ComputeLCM(ComputeLCM(num_process_rows, &
    !      & num_process_columns),num_process_slices),block_multiplier)

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
  PURE SUBROUTINE SplitToLocalBlocks(this, matrix_to_split)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    TYPE(SparseMatrix_t), INTENT(in) :: matrix_to_split
    !! Local Data
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: column_split
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: row_split
    TYPE(SparseMatrix_t) :: tempmat
    INTEGER :: column_counter, row_counter

    !! First Split By Columns
    ALLOCATE(column_split(number_of_blocks_columns))
    ALLOCATE(row_split(number_of_blocks_rows))
    CALL SplitSparseMatrixColumns(matrix_to_split, &
         & number_of_blocks_columns, column_split)

    !! Now Split By Rows
    DO column_counter=1,number_of_blocks_columns
       CALL TransposeSparseMatrix(column_split(column_counter),tempmat)
       CALL SplitSparseMatrixColumns(tempmat, number_of_blocks_rows, &
            & row_split)
       !! And put back in the right place
       DO row_counter=1,number_of_blocks_rows
          CALL TransposeSparseMatrix(row_split(row_counter), &
               & this%local_data(column_counter,row_counter))
       END DO
    END DO

    CALL DestructSparseMatrix(tempmat)
    DO column_counter=1,number_of_blocks_columns
       CALL DestructSparseMatrix(column_split(column_counter))
    END DO
    DO row_counter=1,number_of_blocks_rows
       CALL DestructSparseMatrix(row_split(row_counter))
    END DO
    DEALLOCATE(row_split)
    DEALLOCATE(column_split)
  END SUBROUTINE SplitToLocalBlocks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[in] matrix_to_split the matrix to split up.
  PURE SUBROUTINE MergeLocalBlocks(this, merged_matrix)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    TYPE(SparseMatrix_t), INTENT(inout) :: merged_matrix
    !! Local Data
    TYPE(SparseMatrix_t), DIMENSION(:), ALLOCATABLE :: rows_of_merged
    TYPE(SparseMatrix_t) :: tempmat
    INTEGER :: row_counter

    ALLOCATE(rows_of_merged(number_of_blocks_rows))

    DO row_counter = 1, number_of_blocks_rows
       CALL ComposeSparseMatrixColumns(this%local_data(:,row_counter), tempmat)
       CALL TransposeSparseMatrix(tempmat,rows_of_merged(row_counter))
    END DO

    CALL ComposeSparseMatrixColumns(rows_of_merged,tempmat)
    CALL TransposeSparseMatrix(tempmat,merged_matrix)

    !! Cleanup
    CALL DestructSparseMatrix(tempmat)
    DO row_counter=1,number_of_blocks_rows
       CALL DestructSparseMatrix(rows_of_merged(row_counter))
    END DO
    DEALLOCATE(rows_of_merged)
  END SUBROUTINE MergeLocalBlocks
END MODULE DistributedBlockedSparseMatrixModule
