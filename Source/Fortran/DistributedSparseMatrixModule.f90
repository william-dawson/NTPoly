!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Operations.
MODULE DistributedSparseMatrixModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE DistributedMatrixMemoryPoolModule, ONLY : &
    & DistributedMatrixMemoryPool_t, ConstructDistributedMatrixMemoryPool, &
    & CheckDistributedMemoryPoolValidity
  USE MatrixGatherModule, ONLY : GatherHelper_t, &
       & GatherAndComposeData, GatherANdComposeCleanup, &
       & GatherAndSumData, GatherAndSumCleanup, &
       & GatherSizes
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation
  USE ProcessGridModule, ONLY : grid_error, IsRoot, RootID, &
       & num_process_slices, num_process_rows, num_process_columns, &
       & slice_size, my_slice, my_row, my_column, &
       & within_slice_rank, between_slice_rank, row_rank, column_rank, &
       & global_comm, within_slice_comm, row_comm, column_comm, &
       & between_slice_comm
  USE SparseMatrixModule, ONLY : SparseMatrix_t, ConstructEmptySparseMatrix, &
       & ConstructFromTripletList, TransposeSparseMatrix, SparseMatrixNorm, &
       & ScaleSparseMatrix, PrintSparseMatrix, MatrixToTripletList, &
       & IncrementSparseMatrix, Gemm, DestructSparseMatrix
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletModule, ONLY : Triplet_t, GetMPITripletType
  USE TripletListModule, ONLY : TripletList_t, ConstructTripletList, &
      & DestructTripletList, SortTripletList
  USE mpi
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for a distributed CSR matrix.
  !! A CSR matrix stores locally a block of SparseMatrix types.
  !! That local CSR matrix is also divided up into blocks for sending/receiving
  !! purposes. This allows us to deal with non square slices of a matrix, and
  !! is also important for the 3D multiply algorithm.
  TYPE, PUBLIC :: DistributedSparseMatrix
     !> Number of matrix rows/columns for the full matrix, scaled for process
     !! grid.
     INTEGER :: logical_matrix_dimension
     !> Number of matrix rows/columns for the full matrix, unscaled.
     INTEGER :: actual_matrix_dimension
     !! Local Storage
     TYPE(SparseMatrix_t) :: local_data !< local CSR matrix.
     INTEGER :: start_column !< first column stored locally.
     INTEGER :: end_column !< last column stored locally  is less than this.
     INTEGER :: start_row !< first row stored locally.
     INTEGER :: end_row !< last row stored locally is less than this.
     INTEGER :: local_columns !< number of local columns.
     INTEGER :: local_rows !< number of local rows.
     !! Blocking Information
     INTEGER :: number_of_blocks_columns !< number of column blocks.
     INTEGER :: number_of_blocks_rows !< number of row blocks.
     INTEGER :: size_of_blocks_columns !< number of columns in a block.
     INTEGER :: size_of_blocks_rows !< number of rows in a block.
     !> Who gets which block.
     INTEGER, DIMENSION(:), ALLOCATABLE :: send_start_array_columns
     !> Who gets which block.
     INTEGER, DIMENSION(:), ALLOCATABLE :: send_start_array_rows
     !> Load Balancing Information
     TYPE(Permutation_t) :: BalancePermutation
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
  !! Basic Linear Algebra
  PUBLIC :: IncrementDistributedSparseMatrix
  PUBLIC :: DistributedGemm
  PUBLIC :: ScaleDistributedSparseMatrix
  PUBLIC :: DistributedSparseNorm
  PUBLIC :: Trace
  PUBLIC :: EigenCircle
  PUBLIC :: ComputeSigma
  !! Utilities
  PUBLIC :: PrintDistributedSparseMatrix
  PUBLIC :: GetSize
  PUBLIC :: GetLoadBalance
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix.
  !! @param[out] this the matrix to be constructed.
  !! @param[in] matrix_dim_ the dimension of the full matrix.
  SUBROUTINE ConstructEmpty(this, matrix_dim_)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    INTEGER, INTENT(in) :: matrix_dim_
    !! Some Temporary Variables
    INTEGER :: space_between_blocks_columns, space_between_blocks_rows
    INTEGER :: counter
    !  integer :: alloc_error

    CALL DestructDistributedSparseMatrix(this)

    this%actual_matrix_dimension = matrix_dim_
    this%logical_matrix_dimension = CalculateScaledDimension(matrix_dim_)
    !! Local Block Size Description
    this%local_rows = this%logical_matrix_dimension/num_process_rows
    this%local_columns = this%logical_matrix_dimension/num_process_columns
    !! Which Block Does This Process Hold?
    this%start_row = this%local_rows * my_row + 1
    this%end_row   = this%start_row + this%local_rows
    this%start_column = this%local_columns * my_column + 1
    this%end_column   = this%start_column + this%local_columns

    !! Last Things Is To Build Virtual Blocks For Handling 3D Multiply Case
    !! Columns
    this%number_of_blocks_columns = (num_process_rows/num_process_columns)* &
         & num_process_slices
    IF (this%number_of_blocks_columns .EQ. 0) THEN
       this%number_of_blocks_columns = 1*num_process_slices
    END IF
    this%size_of_blocks_columns = &
         & this%local_columns/this%number_of_blocks_columns
    space_between_blocks_columns = &
         & this%size_of_blocks_columns*num_process_slices
    ALLOCATE(this%send_start_array_columns &
         & (this%number_of_blocks_columns/num_process_slices))
    this%send_start_array_columns(1) = my_slice * this%size_of_blocks_columns +1
    DO counter = 2, this%number_of_blocks_columns/num_process_slices
       this%send_start_array_columns(counter) = &
            & this%send_start_array_columns(counter-1) + &
            & space_between_blocks_columns
    END DO
    !! Rows
    this%number_of_blocks_rows = (num_process_columns/num_process_rows)* &
         & num_process_slices
    IF (this%number_of_blocks_rows .EQ. 0) THEN
       this%number_of_blocks_rows = 1*num_process_slices
    END IF
    this%size_of_blocks_rows = this%local_rows/this%number_of_blocks_rows
    space_between_blocks_rows = this%size_of_blocks_rows*num_process_slices
    ALLOCATE(this%send_start_array_rows &
         & (this%number_of_blocks_rows/num_process_slices))
    this%send_start_array_rows(1) = my_slice * this%size_of_blocks_rows + 1
    DO counter = 2, this%number_of_blocks_rows/num_process_slices
       this%send_start_array_rows(counter) = &
            & this%send_start_array_rows(counter-1) + space_between_blocks_rows
    END DO

    CALL ConstructDefaultPermutation(this%BalancePermutation, &
         & this%logical_matrix_dimension)
  END SUBROUTINE ConstructEmpty
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix
  !! @param[in,out] this the matrix to destruct
  SUBROUTINE DestructDistributedSparseMatrix(this)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this

    CALL DestructSparseMatrix(this%local_data)
    IF (ALLOCATED(this%send_start_array_columns)) THEN
       DEALLOCATE(this%send_start_array_columns)
    END IF
    IF (ALLOCATED(this%send_start_array_rows)) THEN
       DEALLOCATE(this%send_start_array_rows)
    END IF
  END SUBROUTINE DestructDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  SUBROUTINE CopyDistributedSparseMatrix(matA, matB)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in)    :: matA
    TYPE(DistributedSparseMatrix), INTENT(inout) :: matB

    CALL DestructDistributedSparseMatrix(matB)
    matB = matA
  END SUBROUTINE CopyDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !! should pass in triplet lists with global coordinates. It doesn't matter
  !! where each triplet is stored, as long as global coordinates are given.
  !! @param[inout] this the matrix to fill.
  !! @param[in] triplet_list the triplet list of values.
  SUBROUTINE FillFromTripletList(this,triplet_list)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: this
    TYPE(TripletList_t), INTENT(in) :: triplet_list
    !! Local Data
    TYPE(Permutation_t) :: basic_permutation
    TYPE(TripletList_t) :: sorted_triplet_list

    !! First we redistribute the triplet list to get all the local data
    !! on the correct process.
    CALL ConstructDefaultPermutation(basic_permutation, &
         & this%logical_matrix_dimension)
    CALL RedistributeTripletList(this,basic_permutation%index_lookup, &
         & basic_permutation%reverse_index_lookup, triplet_list, &
         & sorted_triplet_list)

    !! Now we can just construct a local matrix.
    CALL ConstructFromTripletList(this%local_data,sorted_triplet_list, &
         & this%local_rows, this%local_columns)

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
    CALL ConstructFromTripletList(this%local_data,sorted_triplet_list, &
         & this%local_rows,this%local_columns)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(unsorted_triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE FillDistributedIdentity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
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
    CALL ConstructFromTripletList(this%local_data,sorted_triplet_list, &
         & this%local_rows,this%local_columns)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(unsorted_triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE FillDistributedPermutation
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
    !! Reading The File
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: unsorted_triplet_list
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
    INTEGER :: triplet_counter

    !! Setup Involves Just The Root Opening And Reading Parameter Data
    CALL StartTimer("MPI Read Text")
    bytes_per_character = sizeof(temp_char)
    IF (IsRoot()) THEN
       header_length = 0
       local_file_handler = 16
       OPEN(local_file_handler, file=file_name, status="old")
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

    !! Build Local Storage
    CALL ConstructEmpty(this,matrix_rows)

    !! Global read
    CALL MPI_File_open(within_slice_comm,file_name,MPI_MODE_RDONLY, &
         & MPI_INFO_NULL,mpi_file_handler,grid_error)

    !! Compute Offsets
    CALL MPI_File_get_size(mpi_file_handler,total_file_size,grid_error)
    local_data_size = (total_file_size - bytes_per_character*header_length)/&
         slice_size
    local_offset = bytes_per_character*header_length + &
         local_data_size*within_slice_rank
    local_data_size_plus_buffer = local_data_size
    IF (.NOT. within_slice_rank .EQ. slice_size-1) THEN
       local_data_size_plus_buffer = local_data_size + &
            MAX_LINE_LENGTH*bytes_per_character
    ELSE
       local_data_size_plus_buffer = (total_file_size - local_offset)
    END IF
    ALLOCATE(CHARACTER(LEN=local_data_size_plus_buffer) :: mpi_input_buffer)

    !! Do Actual Reading
    CALL MPI_File_read_at_all(mpi_file_handler,local_offset,mpi_input_buffer, &
         INT(local_data_size_plus_buffer),MPI_CHARACTER,mpi_status,grid_error)

    !! Trim Off The Half Read Line At The Start
    IF (.NOT. within_slice_rank .EQ. RootID) THEN
       full_buffer_counter = INDEX(mpi_input_buffer,new_line('A')) + 1
    ELSE
       full_buffer_counter = 1
    END IF

    !! Read By Line
    end_of_buffer = .FALSE.
    triplet_counter = 1
    !! WDDAWSON: FIX THIS. We should be able to use less memory.
    CALL ConstructTripletList(triplet_list,this%local_rows*this%local_columns*2)
    DO WHILE(.NOT. end_of_buffer)
       current_line_length = INDEX(mpi_input_buffer(full_buffer_counter:),&
            new_line('A'))

       IF (current_line_length .EQ. 0) THEN !! Hit The End Of The Buffer
          end_of_buffer = .TRUE.
       ELSE
          temp_substring = mpi_input_buffer(full_buffer_counter:&
               & full_buffer_counter+current_line_length-1)
          READ(temp_substring(:current_line_length-1),*) temp_triplet%index_row, &
               & temp_triplet%index_column, temp_triplet%point_value
          triplet_list%data(triplet_counter) = temp_triplet
          triplet_counter = triplet_counter + 1
          IF (full_buffer_counter + current_line_length .GE. &
               & local_data_size+2) THEN
             IF (.NOT. within_slice_rank .EQ. slice_size-1) THEN
                end_of_buffer = .TRUE.
             END IF
          END IF
          full_buffer_counter = full_buffer_counter + current_line_length
       END IF
    END DO

    !! Cleanup
    CALL MPI_File_close(mpi_file_handler,grid_error)
    CALL StopTimer("MPI Read Text")

    !! Redistribute The Matrix
    CALL ConstructTripletList(unsorted_triplet_list,triplet_counter-1)
    unsorted_triplet_list%data = triplet_list%data(:triplet_counter-1)
    CALL FillFromTripletList(this,unsorted_triplet_list)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(unsorted_triplet_list)
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
    INTEGER :: local_file_handler
    INTEGER :: mpi_file_handler
    !! Reading The File
    TYPE(TripletList_t) :: triplet_list
    INTEGER :: matrix_rows, matrix_columns, total_values
    INTEGER :: local_triplets
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER :: bytes_per_int, bytes_per_double
    INTEGER :: triplet_mpi_type
    !! Temporary variables
    INTEGER :: mpi_status(MPI_STATUS_SIZE)

    CALL StartTimer("MPI Read Binary")
    !! Get The Matrix Parameters
    CALL MPI_Type_extent(MPI_INT,bytes_per_int,grid_error)
    CALL MPI_Type_extent(MPINTREAL,bytes_per_double,grid_error)
    IF (IsRoot()) THEN
       local_file_handler = 16
       OPEN(local_file_handler, file=file_name, form="unformatted", &
            access="direct",status="old",recl=bytes_per_int)
       READ(local_file_handler,rec=1) matrix_rows
       READ(local_file_handler,rec=2) matrix_columns
       READ(local_file_handler,rec=3) total_values
       CLOSE(local_file_handler)
    END IF

    !! Broadcast Parameters
    CALL MPI_Bcast(matrix_rows,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(matrix_columns,1,MPI_INT,RootID,global_comm,grid_error)
    CALL MPI_Bcast(total_values,1,MPI_INT,RootID,global_comm,grid_error)

    !! Build Local Storage
    CALL ConstructEmpty(this, matrix_rows)

    !! Compute Offset
    local_triplets = total_values/slice_size
    local_offset = local_triplets * (within_slice_rank)
    header_size = 3 * bytes_per_int
    IF (within_slice_rank .EQ. slice_size - 1) THEN
       local_triplets = INT(total_values) - INT(local_offset)
    END IF
    local_offset = local_offset*(bytes_per_int*2+bytes_per_double) + header_size
    CALL ConstructTripletList(triplet_list,local_triplets)

    !! Do The Actual Reading
    CALL MPI_File_open(within_slice_comm,file_name,MPI_MODE_RDONLY,MPI_INFO_NULL,&
         & mpi_file_handler,grid_error)

    triplet_mpi_type = GetMPITripletType()
    CALL MPI_File_set_view(mpi_file_handler,local_offset,triplet_mpi_type,&
         & triplet_mpi_type,"native",MPI_INFO_NULL,grid_error)
    CALL MPI_File_read_all(mpi_file_handler,triplet_list%data,local_triplets,&
         & triplet_mpi_type,mpi_status,grid_error)
    CALL MPI_File_close(mpi_file_handler,grid_error)
    CALL StopTimer("MPI Read Binary")

    !! Redistribute
    CALL FillFromTripletList(this,triplet_list)
    ! call RedistributeTripletList(this,this%BalancePermutation%index_lookup, &
    !      & this%BalancePermutation%reverse_index_lookup, &
    !      & triplet_list, sorted_triplet_list)
    ! call ConstructFromTripletList(this%local_data,sorted_triplet_list, &
    !      & this%local_rows,this%local_columns)
    CALL DestructTripletList(triplet_list)
    ! call DestructTripletList(sorted_triplet_list)
  END SUBROUTINE ConstructFromBinary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !> Construct a distributed sparse matrix from a dense matrix binary file in
  ! !! parallel.
  ! !! @param[out] this the file being constructed.
  ! !! @param[in] file_name name of the file to read.
  ! !! @param[in] threshold for flushing small matrix values.
  ! !! @param[in] matrix_rows the number of rows in this matrix.
  ! !! @param[in] matrix_columns the number of columns in this matrix.
  ! subroutine ConstructFromDenseBinary(this, file_name, threshold, matrix_rows, &
  !   & matrix_columns)
  !   !! Parameters
  !   type(DistributedSparseMatrix) :: this
  !   character(len=*), intent(in) :: file_name
  !   real(NTREAL), intent(in) :: threshold
  !   integer :: matrix_rows, matrix_columns
  !   !! File Handles
  !   integer :: mpi_file_handler
  !   !! Reading The File
  !   real(NTREAL), dimension(:), allocatable :: double_list
  !   type(Triplet_t), dimension(:), allocatable :: triplet_list
  !   type(Triplet_t), dimension(:), allocatable :: sorted_triplet_list
  !   integer :: local_values
  !   integer :: values_offset
  !   integer(KIND=MPI_OFFSET_KIND) :: local_offset
  !   integer :: bytes_per_double
  !   !! Temporary variables
  !   integer :: counter
  !   integer :: triplet_list_counter
  !   integer :: temp_row, temp_column
  !   integer :: mpi_status(MPI_STATUS_SIZE)
  !
  !   !! Build Local Storage
  !   call ConstructEmpty(this,matrix_rows)
  !
  !   !! Size Information
  !   call MPI_Type_extent(MPINTREAL,bytes_per_double,grid_error)
  !
  !   !! Compute Offset
  !   local_values = this%local_columns * this%local_rows
  !   values_offset = local_values * (within_slice_rank)
  !   local_offset =  bytes_per_double * values_offset
  !   allocate(triplet_list(local_values))
  !   allocate(double_list(local_values))
  !
  !   !! Do Actual Reading
  !   call MPI_File_open(within_slice_comm, file_name, MPI_MODE_RDONLY,   &
  !        & MPI_INFO_NULL, mpi_file_handler, grid_error)
  !   call MPI_File_set_view(mpi_file_handler, local_offset, MPINTREAL, &
  !        & MPINTREAL,"native", MPI_INFO_NULL, grid_error)
  !   call MPI_File_read_all(mpi_file_handler, double_list, local_values, &
  !        & MPINTREAL, mpi_status, grid_error)
  !   call MPI_File_close(mpi_file_handler,grid_error)
  !
  !   !! Create a triplet list
  !   triplet_list_counter = 0
  !   prune_list: do counter = 1, local_values
  !     if (abs(double_list(counter)) .gt. threshold) then
  !       triplet_list_counter = triplet_list_counter + 1
  !       temp_column = mod((values_offset+counter-1),matrix_columns) + 1
  !       temp_row = (values_offset+counter-1)/matrix_columns + 1
  !       triplet_list(triplet_list_counter)%index_column = temp_column
  !       triplet_list(triplet_list_counter)%index_row = temp_row
  !       triplet_list(triplet_list_counter)%point_value = double_list(counter)
  !     end if
  !   end do prune_list
  !
  !   !! Redistribute
  !   call RedistributeTripletList(this,this%BalancePermutation%index_lookup, &
  !        & this%BalancePermutation%reverse_index_lookup, &
  !        & triplet_list(:triplet_list_counter), sorted_triplet_list)
  !   call ConstructFromTripletList(this%local_data,sorted_triplet_list, &
  !        & this%local_rows,this%local_columns)
  !
  !   !! Cleanup
  !   deallocate(double_list)
  !   deallocate(triplet_list)
  !   deallocate(sorted_triplet_list)
  ! end subroutine ConstructFromDenseBinary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !subroutine ConstructDistributedSparseMatrixFromFileSet(this, file_name)
  !  !! Parameters
  !  type(DistributedSparseMatrix) :: this
  !  character(len=*), intent(in) :: file_name
  !  !! Reading The File
  !  character(len=100) :: split_file_name
  !  integer :: split_file_number
  !  character(len=20) :: split_file_number_str
  !  type(Triplet_t), dimension(:), allocatable :: triplet_list
  !  type(Triplet_t), dimension(:), allocatable :: sorted_triplet_list
  !  integer :: file_handler
  !  integer :: matrix_columns, matrix_rows, total_values
  !  integer :: matrix_dim
  !
  !  !! Determine the split file name
  !  split_file_number = my_row * num_process_columns + my_column
  !  write(split_file_number_str,*) split_file_number
  !  split_file_name = trim(file_name)//"-"//&
  !  &   trim(adjustl(split_file_number_str))//".mtx"
  !
  !  file_handler = 16
  !  open(file_handler,file=split_file_name,status='old')
  !  !! Comment Lines
  !  read(file_handler,*)
  !  read(file_handler,*)
  !  !! About The Matrix
  !  read(file_handler,*) matrix_rows, matrix_columns, total_values
  !  call MPI_Allreduce(matrix_rows,matrix_dim,1,MPI_INT,MPI_SUM,row_comm,&
  !       & grid_error)
  !
  !  !! Build Local Storage
  !  call ConstructEmptySparseDistributedMatrix(this,matrix_rows, matrix_columns,&
  !       & matrix_dim)
  !
  !  !! Finish Reading
  !  allocate(triplet_list(this%local_rows*this%local_columns))
  !  call ReadInTripletList(file_handler, total_values, triplet_list)
  !  call SortTripletList(triplet_list(:total_values),this%local_columns,&
  !       & sorted_triplet_list)
  !  call ConstructFromTripletList(this%local_data,sorted_triplet_list, &
  !       & this%local_rows,this%local_columns)
  !
  !  !! Finished With File
  !  close(file_handler)
  !  deallocate(triplet_list)
  !  deallocate(sorted_triplet_list)
  !end subroutine ConstructDistributedSparseMatrixFromFileSet
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

    !! Determine Write Location
    bytes_per_int = sizeof(temp_int)
    bytes_per_double = sizeof(temp_double)
    local_data_size = SIZE(this%local_data%values)*(bytes_per_int*2 + &
         & bytes_per_double*1)
    header_size = bytes_per_int*3
    ALLOCATE(local_values_buffer(slice_size))
    CALL MPI_Allgather(SIZE(this%local_data%values),1,MPI_INT,&
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

       CALL MatrixToTripletList(this%local_data, triplet_list)
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
            & SIZE(triplet_list%data),triplet_mpi_type,MPI_STATUS_IGNORE, &
            & grid_error)

       !! Cleanup
       CALL MPI_File_close(mpi_file_handler,grid_error)
       CALL DestructTripletList(triplet_list)
    END IF
    DEALLOCATE(local_values_buffer)
    CALL MPI_Barrier(global_comm,grid_error)
  END SUBROUTINE WriteToBinary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a matrix market file.
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
    TYPE(TripletList_t)  :: triplet_list
    INTEGER :: triplet_list_string_length
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER(KIND=MPI_OFFSET_KIND) :: write_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_offset
    INTEGER(KIND=MPI_OFFSET_KIND), PARAMETER :: zero_size = 0
    !! Strings
    CHARACTER(len=*), PARAMETER :: header_line1 = &
         & "%%MatrixMarket matrix coordinate real general" &
         & //new_line('A')//"%"//new_line('A')
    CHARACTER(len=:), ALLOCATABLE :: header_line2
    CHARACTER(len=:), ALLOCATABLE :: write_buffer
    !! Temporary Values
    INTEGER :: counter
    INTEGER :: offset_counter
    INTEGER, PARAMETER :: NEW_LINE_LENGTH = LEN(new_line('A'))
    CHARACTER(len=MAX_LINE_LENGTH) :: temp_string
    INTEGER :: temp_length
    INTEGER :: bytes_per_character
    CHARACTER(len=1) :: temp_char

    bytes_per_character = sizeof(temp_char)

    !! Create the matrix size line
    WRITE(temp_string,*) this%actual_matrix_dimension, &
         & this%actual_matrix_dimension, GetSize(this)
    ALLOCATE(CHARACTER(len=LEN_TRIM(temp_string)+NEW_LINE_LENGTH) :: header_line2)
    WRITE(header_line2,*) TRIM(temp_string)
    header_line2 = header_line2//new_line('A')

    header_size = LEN(header_line1) + LEN(header_line2)

    !! Local Data
    CALL MatrixToTripletList(this%local_data, triplet_list)

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
       WRITE(temp_string,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value, &
            & new_line('A')
       triplet_list_string_length = triplet_list_string_length + &
            & LEN_TRIM(temp_string)
       triplet_list_string_length = triplet_list_string_length + NEW_LINE_LENGTH
    END DO

    !! Write that string to the write buffer
    ALLOCATE(CHARACTER(len=triplet_list_string_length+1) :: write_buffer)
    offset_counter = 1
    DO counter = 1, triplet_list%CurrentSize
       WRITE(temp_string,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value, &
            & new_line('A')
       temp_length = LEN_TRIM(temp_string)+NEW_LINE_LENGTH
       WRITE(write_buffer(offset_counter:),*) temp_string(1:temp_length)
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
            & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL,mpi_file_handler,&
            & grid_error)
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
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !! This will utilize the sparse vector increment routine.
  !! @param[in] matA Matrix A
  !! @param[in,out] matB Matrix B
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

    CALL IncrementSparseMatrix(matA%local_data, matB%local_data, alpha, &
         & threshold)
  END SUBROUTINE IncrementDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  !! @param[inout] matA Matrix A.
  !! @param[in] constant scale factor.
  PURE SUBROUTINE ScaleDistributedSparseMatrix(matA,constant)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(inout) :: matA
    REAL(NTREAL), INTENT(in) :: constant

    CALL ScaleSparseMatrix(matA%local_data,constant)
  END SUBROUTINE ScaleDistributedSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !! C := alpha*matA*matB+ beta*matC
  !! @param[in] matA Matrix A
  !! @param[in] matB Matrix B
  !! @param[in] alpha_in scales the multiplication
  !! @param[in] beta_in scales matrix we sum on to
  !! @param[out] matC = alpha*matA*matB + beta*matC
  !! @param[in] threshold_in for flushing values to zero. Default value is 0.0.
  !! @param[inout] memory_pool_in a memory pool that can be used for the
  !! calculation.
  SUBROUTINE DistributedGemm(matA,matB,matC,alpha_in,beta_in,threshold_in,&
       & memory_pool_in)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: matA
    TYPE(DistributedSparseMatrix), INTENT(in) :: matB
    REAL(NTREAL), OPTIONAL, INTENT(in) :: alpha_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: beta_in
    REAL(NTREAL), OPTIONAL, INTENT(in) :: threshold_in
    TYPE(DistributedMatrixMemoryPool_t), OPTIONAL, INTENT(inout) :: memory_pool_in
    TYPE(DistributedSparseMatrix), INTENT(inout) :: matC
    !! Local Composed Matrix
    TYPE(SparseMatrix_t) :: local_column_block, local_row_block
    TYPE(SparseMatrix_t) :: local_column_blockT
    TYPE(SparseMatrix_t) :: gathered_column_block, gathered_row_block
    TYPE(SparseMatrix_t) :: multiplied_intermediate
    TYPE(TripletList_t) :: full_triplet_list1
    TYPE(TripletList_t) :: full_triplet_list2
    TYPE(TripletList_t) :: triplet_list_column
    TYPE(TripletList_t) :: triplet_list_row
    TYPE(TripletList_t) :: temporary_list
    !! Communication Helpers
    TYPE(GatherHelper_t) :: row_helper, column_helper, sum_helper
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    !! Temporary Variables
    TYPE(DistributedSparseMatrix) :: matAB
    INTEGER :: counter
    INTEGER :: inner_counter
    INTEGER :: running_counter
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold

    CALL StartTimer("GEMM")
    !! Handle the parameters
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

    CALL DestructSparseMatrix(matC%local_data)

    CALL ConstructEmpty(matAB,matA%actual_matrix_dimension)

    !! Triplet Value Based Method To Perform A Local Gather
    CALL MatrixToTripletList(matA%local_data, full_triplet_list1)
    CALL ConstructTripletList(triplet_list_row,full_triplet_list1%CurrentSize)
    running_counter = 1
    DO counter = 1, full_triplet_list1%CurrentSize
       DO inner_counter = 1, SIZE(matA%send_start_array_columns)
          IF (full_triplet_list1%data(counter)%index_column .GE. &
               & matA%send_start_array_columns(inner_counter) .AND. &
               & full_triplet_list1%data(counter)%index_column .LT. &
               & matA%send_start_array_columns(inner_counter)+ &
               & matA%size_of_blocks_columns) THEN
             triplet_list_row%data(running_counter) = &
                  & full_triplet_list1%data(counter)
             running_counter = running_counter + 1
          END IF
       END DO
    END DO
    CALL DestructTripletList(full_triplet_list1)
    CALL ConstructTripletList(temporary_list,running_counter-1)
    temporary_list%data = triplet_list_row%data(:running_counter-1)
    CALL ConstructFromTripletList(local_row_block, &
         & temporary_list, &
         & matA%local_data%rows, &
         & matA%local_data%columns)
    CALL DestructTripletList(triplet_list_row)
    CALL DestructTripletList(temporary_list)

    CALL MatrixToTripletList(matB%local_data, full_triplet_list2)
    CALL ConstructTripletList(triplet_list_column, &
         & SIZE(full_triplet_list2%data))
    running_counter = 1
    DO counter = 1, full_triplet_list2%CurrentSize
       DO inner_counter = 1, UBOUND(matB%send_start_array_rows,1)
          IF (full_triplet_list2%data(counter)%index_row .GE. &
               & matB%send_start_array_rows(inner_counter) .AND. &
               & full_triplet_list2%data(counter)%index_row .LT. &
               & matB%send_start_array_rows(inner_counter)+&
               & matB%size_of_blocks_rows) THEN
             triplet_list_column%data(running_counter) = &
                  & full_triplet_list2%data(counter)
             running_counter = running_counter + 1
          END IF
       END DO
    END DO
    CALL DestructTripletList(full_triplet_list2)
    CALL ConstructTripletList(temporary_list,running_counter-1)
    temporary_list%data = triplet_list_column%data(:running_counter-1)
    CALL ConstructFromTripletList(local_column_block, &
         & temporary_list, &
         & matB%local_data%rows,&
         & matB%local_data%columns)
    CALL DestructTripletList(triplet_list_column)
    CALL DestructTripletList(temporary_list)

    !! Next Is The Global AllGather
    CALL StartTimer("Gather1")
    CALL GatherSizes(local_row_block,row_comm,row_helper)
    CALL MPI_Wait(row_helper%size_request,mpi_status,grid_error)
    CALL GatherAndComposeData(local_row_block,row_comm,gathered_row_block, &
         & row_helper)
    CALL StopTimer("Gather1")
    CALL TransposeSparseMatrix(local_column_block, local_column_blockT)
    CALL DestructSparseMatrix(local_column_block)
    CALL StartTimer("Gather2")
    CALL GatherSizes(local_column_blockT,column_comm,column_helper)
    CALL MPI_Wait(column_helper%size_request,mpi_status,grid_error)
    CALL GatherAndComposeData(local_column_blockT,column_comm, &
         & gathered_column_block, column_helper)
    CALL StopTimer("Gather2")

    CALL StartTimer("Gather1")
    CALL MPI_Wait(row_helper%outer_request,mpi_status,grid_error)
    CALL MPI_Wait(row_helper%inner_request,mpi_status,grid_error)
    CALL MPI_Wait(row_helper%data_request,mpi_status,grid_error)
    CALL GatherAndComposeCleanup(local_row_block,gathered_row_block,row_helper)
    CALL DestructSparseMatrix(local_row_block)
    CALL StopTimer("Gather1")
    CALL StartTimer("Gather2")
    CALL MPI_Wait(column_helper%outer_request,mpi_status,grid_error)
    CALL MPI_Wait(column_helper%inner_request,mpi_status,grid_error)
    CALL MPI_Wait(column_helper%data_request,mpi_status,grid_error)
    CALL GatherAndComposeCleanup(local_column_blockT, &
         & gathered_column_block, column_helper)
    CALL DestructSparseMatrix(local_column_blockT)
    CALL StopTimer("Gather2")

    !! Local Multiply
    IF (PRESENT(memory_pool_in)) THEN
       IF (.NOT. CheckDistributedMemoryPoolValidity(memory_pool_in)) THEN
          CALL ConstructDistributedMatrixMemoryPool(memory_pool_in)
       END IF
    END IF
    CALL StartTimer("Multiply Block")
    !! Rows for gathered_column_blocks because of the transpose
    CALL ConstructEmptySparseMatrix(multiplied_intermediate,&
         & gathered_column_block%rows, gathered_row_block%rows)
    !! Divide by 10, one extra order of magnitude to make sure we don't
    !! prune anything that might cancel out
    IF (num_process_slices .GT. 1) THEN
       IF (.NOT. PRESENT(memory_pool_in)) THEN
          CALL Gemm(gathered_row_block,gathered_column_block,&
               & multiplied_intermediate,IsBTransposed_in=.TRUE.,alpha_in=alpha,&
               & threshold_in=threshold/(num_process_slices*1000))
       ELSE
          CALL Gemm(gathered_row_block,gathered_column_block,&
               & multiplied_intermediate,IsBTransposed_in=.TRUE.,alpha_in=alpha,&
               & threshold_in=threshold/(num_process_slices*1000), &
               & blocked_memory_pool_in=memory_pool_in%grid(1,1))
       END IF
    ELSE
       IF (.NOT. PRESENT(memory_pool_in)) THEN
          CALL Gemm(gathered_row_block,gathered_column_block,   &
               & multiplied_intermediate,IsBTransposed_in=.TRUE.,alpha_in=alpha, &
               & threshold_in=threshold)
       ELSE
          CALL Gemm(gathered_row_block,gathered_column_block,   &
               & multiplied_intermediate,IsBTransposed_in=.TRUE.,alpha_in=alpha, &
               & threshold_in=threshold, &
               & blocked_memory_pool_in=memory_pool_in%grid(1,1))
       END IF
    END IF
    CALL DestructSparseMatrix(gathered_row_block)
    CALL DestructSparseMatrix(gathered_column_block)
    CALL StopTimer("Multiply Block")

    !! Now Sum Along Slices
    CALL StartTimer("Gather3")
    !!call GatherAllAndSum(multiplied_intermediate,between_slice_comm,&
    !!     & matAB%local_data,threshold)
    CALL GatherSizes(multiplied_intermediate,between_slice_comm,sum_helper)
    CALL MPI_Wait(sum_helper%size_request,mpi_status,grid_error)
    CALL GatherAndSumData(multiplied_intermediate,between_slice_comm,sum_helper)
    CALL MPI_Wait(sum_helper%outer_request,mpi_status,grid_error)
    CALL MPI_Wait(sum_helper%inner_request,mpi_status,grid_error)
    CALL MPI_Wait(sum_helper%data_request,mpi_status,grid_error)
    CALL GatherAndSumCleanup(multiplied_intermediate,matAB%local_data, &
         & threshold, sum_helper)
    CALL StopTimer("Gather3")

    CALL CopyDistributedSparseMatrix(matAB,matC)

    CALL DestructDistributedSparseMatrix(matAB)

    CALL DestructSparseMatrix(multiplied_intermediate)
    CALL StopTimer("GEMM")

  END SUBROUTINE DistributedGemm
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

    !! Sum Along Columns
    CALL SparseMatrixNorm(this%local_data,local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
         & MPINTREAL, MPI_SUM,column_comm,grid_error)

    !! Find Max Value Amonst Columns
    norm_value = MAXVAL(local_norm)
    CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
         & row_comm, grid_error)

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

    !! Compute The Local Contribution
    trace_value = 0
    CALL MatrixToTripletList(this%local_data,triplet_list)
    DO counter = 1, triplet_list%CurrentSize
       IF (this%start_row + triplet_list%data(counter)%index_row .EQ. &
            & this%start_column + triplet_list%data(counter)%index_column) THEN
          trace_value = trace_value + triplet_list%data(counter)%point_value
       END IF
    END DO

    !! Sum Among Process Slice
    CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
         & MPI_SUM, within_slice_comm,grid_error)

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

    !! Allocate Space For Result
    ALLOCATE(per_column_min(this%local_data%columns))
    ALLOCATE(per_column_max(this%local_data%columns))

    !! Compute The Local Contribution
    per_column_min = 0
    per_column_max = 0
    CALL MatrixToTripletList(this%local_data,triplet_list)
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
  END SUBROUTINE EigenCircle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print ouf a distributed sparse matrix.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintDistributedSparseMatrix(this, file_name_in)
    !! Parameters
    TYPE(DistributedSparseMatrix) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: file_name_in
    !! The Root Process In Each Column Will Fill These
    TYPE(SparseMatrix_t) :: local_transpose
    TYPE(SparseMatrix_t) :: gathered_matrix_column
    !! The Root Process Will Fill These
    TYPE(SparseMatrix_t) :: gathered_matrix_columnT
    TYPE(SparseMatrix_t) :: gathered_matrix_row

    !! First gather along the columns, which requires first transposing the
    !! local data
    CALL TransposeSparseMatrix(this%local_data, local_transpose)
    CALL GatherAndCompose(local_transpose, column_comm, gathered_matrix_column)

    IF (column_rank .EQ. RootID) THEN
       !! Then gather along the row
       CALL TransposeSparseMatrix(gathered_matrix_column, gathered_matrix_columnT)
       CALL GatherAndCompose(gathered_matrix_columnT, row_comm,gathered_matrix_row)
       !! Make these changes so that it prints the logical rows/columns
       gathered_matrix_row%rows = this%actual_matrix_dimension
       gathered_matrix_row%columns = this%actual_matrix_dimension
       IF (row_rank .EQ. RootID .AND. between_slice_rank .EQ. RootID) THEN
          IF (PRESENT(file_name_in)) THEN
             CALL PrintSparseMatrix(gathered_matrix_row, file_name_in)
          ELSE
             CALL PrintSparseMatrix(gathered_matrix_row)
          END IF
       END IF
    END IF

  END SUBROUTINE PrintDistributedSparseMatrix
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
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter

    ALLOCATE(column_sigma_contribution(this%local_data%columns))
    column_sigma_contribution = 0
    DO outer_counter = 1, this%local_data%columns
       DO inner_counter = this%local_data%outer_index(outer_counter), &
            & this%local_data%outer_index(outer_counter+1)-1
          column_sigma_contribution(outer_counter) = &
               & column_sigma_contribution(outer_counter) + &
               & ABS(this%local_data%values(inner_counter+1))
       END DO
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE,column_sigma_contribution,&
         & this%local_data%columns,MPINTREAL,MPI_SUM, &
         & column_comm, grid_error)
    CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
         & MPINTREAL,MPI_MAX,row_comm,grid_error)
    sigma_value = 1.0d+0/(sigma_value**2)

    DEALLOCATE(column_sigma_contribution)
  END SUBROUTINE ComputeSigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gather all the matrices on to the root process and compose to a big matrix
  !! This is different from GatherAllAndCompose in that only the root process
  !! will receive the data.
  !! @param[in] matrix the local matrix to contribute.
  !! @param[in] communicator mpi communicator along which to gather. Note that
  !! I've marked this in because that's the real semantics, but MPI requires
  !! that it be inout.
  !! @param[out] gathered_matrix = Gathered up matrix. This is not empty only on
  !! the root process.
  SUBROUTINE GatherAndCompose(matrix, communicator, gathered_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(in) :: matrix
    INTEGER, INTENT(inout) :: communicator
    TYPE(SparseMatrix_t) :: gathered_matrix
    !! MPI Information
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_process
    INTEGER, DIMENSION(:), ALLOCATABLE :: displacement
    INTEGER :: counter, inner_counter
    INTEGER :: comm_size, comm_rank
    INTEGER :: temp_offset
    INTEGER :: temporary_total_values

    !! Gather Information About Other Processes
    CALL MPI_Comm_rank(communicator,comm_rank,grid_error)
    CALL MPI_Comm_size(communicator,comm_size,grid_error)
    ALLOCATE(values_per_process(comm_size))
    ALLOCATE(displacement(comm_size))
    CALL MPI_Gather(SIZE(matrix%values),1,MPI_INT,&
         & values_per_process,1,MPI_INT,RootID,communicator,grid_error)
    displacement(1) = 0
    DO counter = 2, SIZE(displacement)
       displacement(counter) = displacement(counter-1) +&
            & values_per_process(counter-1)
    END DO

    !! Allocate Storage for Root Process
    IF (comm_rank == RootID) THEN
       CALL ConstructEmptySparseMatrix(gathered_matrix, &
            & matrix%columns*comm_size,matrix%rows)
       temporary_total_values = SUM(values_per_process)
       ALLOCATE(gathered_matrix%values(temporary_total_values))
       ALLOCATE(gathered_matrix%inner_index(temporary_total_values))
       gathered_matrix%outer_index(1) = 0
    ELSE
       !! This is silly, but we need to allocate the outer index so that
       !! the gather statement (2:) part doesn't crash.
       ALLOCATE(gathered_matrix%outer_index(2))
    END IF

    !! Gather Everything
    CALL MPI_Gatherv(matrix%values,SIZE(matrix%values),MPINTREAL, &
         & gathered_matrix%values, values_per_process, displacement, &
         & MPINTREAL, RootID, communicator, grid_error)
    CALL MPI_Gatherv(matrix%inner_index,SIZE(matrix%values),MPI_INT, &
         & gathered_matrix%inner_index, values_per_process, displacement, &
         & MPI_INT, RootID, communicator, grid_error)
    CALL MPI_Gather(matrix%outer_index(2:), matrix%columns,&
         & MPI_INT, gathered_matrix%outer_index(2:), &
         & matrix%columns, MPI_INT, RootID, communicator, grid_error)

    !! Sum Up The Outer Indices
    IF (comm_rank == RootID) THEN
       DO counter = 1, comm_size - 1
          temp_offset = counter*matrix%columns+1
          DO inner_counter = 1, matrix%columns
             gathered_matrix%outer_index(temp_offset+inner_counter) = &
                  & gathered_matrix%outer_index(temp_offset) + &
                  & gathered_matrix%outer_index(temp_offset+inner_counter)
          END DO
       END DO
    END IF
    DEALLOCATE(values_per_process)
    DEALLOCATE(displacement)

  END SUBROUTINE GatherAndCompose
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
  SUBROUTINE RedistributeTripletList(this,index_lookup,reverse_index_lookup, &
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
       location_list_within_slice(counter)= &
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
  !> Get the total number of non zero entries in the distributed sparse matrix.
  !! @param[in] this the distributed sparse matrix to calculate the non-zero
  !! entries of.
  !! @return the number of non-zero entries in the matrix.
  FUNCTION GetSize(this) RESULT(total_size)
    !! Parameters
    TYPE(DistributedSparseMatrix), INTENT(in) :: this
    INTEGER :: total_size
    !! Local Data
    INTEGER :: local_size

    local_size = SIZE(this%local_data%values)
    CALL MPI_Allreduce(local_size,total_size,1,MPI_INT,MPI_SUM,&
         & within_slice_comm, grid_error)
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

    local_size = SIZE(this%local_data%values)
    CALL MPI_Allreduce(local_size,max_size,1,MPI_INT,MPI_MAX,&
         & within_slice_comm, grid_error)
    CALL MPI_Allreduce(local_size,min_size,1,MPI_INT,MPI_MIN,&
         & within_slice_comm, grid_error)
  END SUBROUTINE GetLoadBalance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate a matrix size that can be divided by the number of processors.
  !! In fact, there are ways to solve this problem which results in smaller
  !! matrices.
  !! @param[in] matrix_dim the dimension of the actual matrix.
  !! @return a new dimension which includes padding.
  !! @todo write a more optimal algorithm using either plain brute force,
  !! or a prime factor based algorithm.
  FUNCTION CalculateScaledDimension(matrix_dim) RESULT(scaled_dim)
    !! Parameters
    INTEGER, INTENT(in) :: matrix_dim
    INTEGER :: scaled_dim
    !! Local Data
    INTEGER :: size_ratio
    INTEGER :: lcm

    lcm = num_process_rows*num_process_columns*num_process_slices

    size_ratio = matrix_dim/lcm
    IF (size_ratio * lcm .EQ. matrix_dim) THEN
       scaled_dim = matrix_dim
    ELSE
       scaled_dim = (size_ratio + 1)*lcm
    END IF
  END FUNCTION CalculateScaledDimension
END MODULE DistributedSparseMatrixModule
