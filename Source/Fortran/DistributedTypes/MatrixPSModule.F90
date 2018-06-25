!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Operations.
MODULE MatrixPSModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE ISO_C_BINDING, ONLY : c_long
  USE LoggingModule
  USE MatrixMarketModule, ONLY : ParseMMHeader
  USE MatrixPModule, ONLY : Matrix_p
  USE MatrixReduceModule, ONLY : ReduceHelper_t
  USE MatrixSModule, ONLY : Matrix_ls, Matrix_lsr, Matrix_lsc, &
      & ComposeMatrix, SplitMatrix
  USE PermutationModule, ONLY : Permutation_t
  USE ProcessGridModule, ONLY : ProcessGrid_t, global_grid
  USE TripletListModule, ONLY : TripletList, TripletList_r, TripletList_c
  USE TripletModule, ONLY : Triplet, Triplet_r, Triplet_c
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, ABSTRACT, EXTENDS(Matrix_p), PUBLIC :: Matrix_ps
  CONTAINS
    !! Fill In Special Matrices
    PROCEDURE(FillMatrixPermutation_ps), DEFERRED :: FillPermutation
    !! Utilities
    PROCEDURE(GetMatrixLoadBalance_ps), DEFERRED :: GetLoadBalance
    PROCEDURE(GetMatrixSize_ps), DEFERRED :: GetSize
    PROCEDURE(FilterMatrix_ps), DEFERRED :: Filter
    PROCEDURE(MergeMatrixLocalBlocks_ps), DEFERRED :: MergeBlocks
    PROCEDURE(SplitMatrixToLocalBlocks_ps), DEFERRED :: SplitToLocal
  END TYPE Matrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE, EXTENDS(Matrix_ps), PUBLIC :: Matrix_psr
     !> A 2D array of local CSR matrices.
     TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: local_data
  CONTAINS
     !! Construct/Destruct
     PROCEDURE :: InitEmpty => ConstructEmptyMatrix_psr
     PROCEDURE :: Destruct => DestructMatrix_psr
     PROCEDURE :: Copy => CopyMatrix_psr
     !! File I/O
     PROCEDURE :: InitMatrixMarket => ConstructFromMatrixMarket_psr
     PROCEDURE :: InitBinary => ConstructFromBinary_psr
     PROCEDURE :: WriteMatrixMarket => WriteMatrixToMatrixMarket_psr
     PROCEDURE :: WriteBinary => WriteMatrixToBinary_psr
     !! Fill In Special Matrices
     PROCEDURE :: FillFromTripletList => FillMatrixFromTripletList_psr
     PROCEDURE :: FillIdentity => FillMatrixIdentity_psr
     PROCEDURE :: FillPermutation => FillMatrixPermutation_psr
     !! Basic Accessor
     PROCEDURE :: GetTripletList => GetTripletList_psr
     PROCEDURE :: GetBlock => GetBlock_psr
     !! Printing To The Console
     PROCEDURE :: Print => PrintMatrix_psr
     PROCEDURE :: PrintInfo => PrintMatrixInformation_psr
     !! Utilities
     PROCEDURE :: GetLoadBalance => GetMatrixLoadBalance_psr
     PROCEDURE :: GetSize => GetMatrixSize_psr
     PROCEDURE :: Filter => FilterMatrix_psr
     PROCEDURE :: MergeBlocks => MergeMatrixLocalBlocks_psr
     PROCEDURE :: SplitToLocal => SplitMatrixToLocalBlocks_psr
     PROCEDURE :: Transpose => TransposeMatrix_psr
     PROCEDURE :: CommSplit => CommSplitMatrix_psr
  END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ABSTRACT INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE FillMatrixPermutation_ps(this, permutation_vector, permuterows)
      IMPORT :: Matrix_ps
      IMPLICIT NONE
      CLASS(Matrix_ps), INTENT(INOUT) :: this
      INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
      LOGICAL, OPTIONAL, INTENT(IN) :: permuterows
   END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE GetMatrixLoadBalance_ps(this, min_size, max_size)
      IMPORT :: Matrix_ps
      IMPLICIT NONE
      CLASS(Matrix_ps), INTENT(IN) :: this
      INTEGER, INTENT(OUT) :: min_size
      INTEGER, INTENT(OUT) :: max_size
   END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   FUNCTION GetMatrixSize_ps(this) RESULT(total_size)
      USE ISO_C_BINDING, ONLY : c_long
      IMPORT :: Matrix_ps
      IMPLICIT NONE
      CLASS(Matrix_ps), INTENT(IN) :: this
      INTEGER(c_long) :: total_size
   END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE FilterMatrix_ps(this, threshold)
     USE DataTypesModule, ONLY : NTREAL
      IMPORT :: Matrix_ps
      IMPLICIT NONE
      CLASS(Matrix_ps), INTENT(INOUT) :: this
      REAL(NTREAL), INTENT(IN) :: threshold
   END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[in] matrix_to_split the matrix to split up.
  SUBROUTINE SplitMatrixToLocalBlocks_ps(this, matrix_to_split)
     USE MatrixSModule, ONLY : Matrix_ls
     IMPORT :: Matrix_ps
     IMPLICIT NONE
     CLASS(Matrix_ps), INTENT(INOUT) :: this
     CLASS(Matrix_ls), INTENT(IN) :: matrix_to_split
  END SUBROUTINE SplitMatrixToLocalBlocks_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE MergeMatrixLocalBlocks_ps(this, merged_matrix)
     USE MatrixSModule, ONLY : Matrix_ls
     IMPORT :: Matrix_ps
     IMPLICIT NONE
     CLASS(Matrix_ps), INTENT(IN) :: this
     CLASS(Matrix_ls), INTENT(INOUT) :: merged_matrix
  END SUBROUTINE MergeMatrixLocalBlocks_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix.
  !! @param[out] this the matrix to be constructed.
  !! @param[in] matrix_dim the dimension of the full matrix.
  !! @param[in] process_grid_in a process grid to host the matrix (optional).
  SUBROUTINE ConstructEmptyMatrix_psr(this, matrix_dim, process_grid_in)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)             :: matrix_dim
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: process_grid_in
    !! Local Variables
    TYPE(Matrix_lsr) :: zeromatrix

    CALL this%Destruct

    !! Process Grid
    IF (PRESENT(process_grid_in)) THEN
       CALL this%process_grid%Copy(process_grid_in)
    ELSE
       CALL this%process_grid%Copy(global_grid)
    END IF

    !! Matrix Dimensions
    this%actual_matrix_dimension = matrix_dim
    this%logical_matrix_dimension = CalculateScaledDimension(this%process_grid,&
        matrix_dim)

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
    ALLOCATE(this%local_data(this%process_grid%number_of_blocks_rows, &
         & this%process_grid%number_of_blocks_columns))
    CALL zeromatrix%InitEmpty(this%local_rows, this%local_columns)

    CALL this%SplitToLocal(zeromatrix)
    CALL zeromatrix%Destruct
  END SUBROUTINE ConstructEmptyMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix.
  !! @param[inout] this the matrix to destruct.
  PURE SUBROUTINE DestructMatrix_psr(this)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: II, JJ

    IF (ALLOCATED(this%local_data)) THEN
       DO JJ = 1, this%process_grid%number_of_blocks_columns
          DO II = 1, this%process_grid%number_of_blocks_rows
             CALL this%local_data(II,JJ)%Destruct
          END DO
       END DO
       DEALLOCATE(this%local_data)
    END IF
  END SUBROUTINE DestructMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] this = matA
  SUBROUTINE CopyMatrix_psr(this, matA)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CLASS(Matrix_p), INTENT(IN)    :: matA
    INTEGER :: II, JJ

    CALL this%Destruct
    SELECT TYPE(MatA)
    CLASS IS (Matrix_psr)
      !! Basic Values
      this%logical_matrix_dimension = matA%logical_matrix_dimension
      this%actual_matrix_dimension = matA%actual_matrix_dimension
      this%start_column = matA%start_column
      this%end_column = matA%end_column
      this%start_row = matA%start_row
      this%end_row = matA%end_row
      this%local_columns = matA%local_columns
      this%local_rows = matA%local_rows
      CALL this%process_grid%copy(matA%process_grid)

      !! Local Data
      ALLOCATE(this%local_data(this%process_grid%number_of_blocks_rows, &
           & this%process_grid%number_of_blocks_columns))
      DO JJ = 1, this%process_grid%number_of_blocks_columns
         DO II = 1, this%process_grid%number_of_blocks_rows
            CALL this%local_data(II,JJ)%Copy(matA%local_data(II,JJ))
         END DO
      END DO

    END SELECT
  END SUBROUTINE CopyMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  !! Read \cite boisvert1996matrix for the details.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructFromMatrixMarket_psr(this, file_name)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
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
    bytes_per_character = sizeof(temp_char)
    IF (global_grid%IsRoot()) THEN
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
    CALL this%InitEmpty(matrix_rows, global_grid)

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
    CALL triplet_list%Init
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
             CALL triplet_list%Append(temp_triplet)
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
    CALL triplet_list%Symmetrize(pattern_type)
    CALL this%FillFromTripletList(triplet_list)

    CALL triplet_list%Destruct
    DEALLOCATE(mpi_input_buffer)
  END SUBROUTINE ConstructFromMatrixMarket_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  !! Faster than text, so this is good for check pointing.
  !! @param[out] this the file being constructed.
  !! @param[in] file_name name of the file to read.
  SUBROUTINE ConstructFromBinary_psr(this, file_name)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    ! !! File Handles
    ! INTEGER :: mpi_file_handler
    ! !! Reading The File
    ! TYPE(TripletList_r) :: triplet_list
    ! INTEGER :: matrix_rows, matrix_columns, total_values
    ! INTEGER, DIMENSION(3) :: matrix_information
    ! INTEGER :: local_triplets
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    ! INTEGER :: bytes_per_int, bytes_per_double
    ! INTEGER :: triplet_mpi_type
    ! !! Temporary variables
    ! INTEGER :: mpi_status(MPI_STATUS_SIZE)
    ! INTEGER :: ierr
    ! TYPE(Triplet_r) :: temp_triplet
    !
    ! TYPE(Triplet_r), DIMENSION(2) :: test
    !
    ! CALL MPI_Type_extent(MPI_INT,bytes_per_int,ierr)
    ! CALL MPI_Type_extent(MPINTREAL,bytes_per_double,ierr)
    !
    ! CALL MPI_File_open(global_grid%global_comm,file_name,MPI_MODE_RDONLY,&
    !      & MPI_INFO_NULL,mpi_file_handler,ierr)
    ! IF (ierr .NE. 0) THEN
    !    IF (global_grid%IsRoot()) THEN
    !       WRITE(*,*) file_name, " doesn't exist"
    !    END IF
    !    CALL MPI_Abort(global_grid%global_comm, -1, ierr)
    ! END IF
    !
    ! !! Get The Matrix Parameters
    ! IF (global_grid%IsRoot()) THEN
    !    local_offset = 0
    !    CALL MPI_File_read_at(mpi_file_handler, local_offset, &
    !         & matrix_information, 3, MPI_INT, mpi_status, ierr)
    !    matrix_rows = matrix_information(1)
    !    matrix_columns = matrix_information(2)
    !    total_values = matrix_information(3)
    ! END IF
    !
    ! !! Broadcast Parameters
    ! CALL MPI_Bcast(matrix_rows, 1, MPI_INT, global_grid%RootID, &
    !      & global_grid%global_comm, ierr)
    ! CALL MPI_Bcast(matrix_columns, 1, MPI_INT, global_grid%RootID, &
    !      & global_grid%global_comm, ierr)
    ! CALL MPI_Bcast(total_values, 1, MPI_INT ,global_grid%RootID, &
    !      & global_grid%global_comm, ierr)
    !
    ! !! Build Local Storage
    ! CALL this%InitEmpty(matrix_rows, global_grid)
    !
    ! !! Compute Offset
    ! local_triplets = total_values/this%process_grid%total_processors
    ! local_offset = local_triplets * (this%process_grid%global_rank)
    ! header_size = 3 * bytes_per_int
    ! IF (this%process_grid%global_rank .EQ. &
    !      & this%process_grid%total_processors - 1) THEN
    !    local_triplets = INT(total_values) - INT(local_offset)
    ! END IF
    ! local_offset = local_offset*(bytes_per_int*2+bytes_per_double) + header_size
    ! CALL triplet_list%Init(local_triplets)
    !
    ! !! Do The Actual Reading
    ! triplet_mpi_type = temp_triplet%GetMPIType()
    ! CALL MPI_File_set_view(mpi_file_handler,local_offset,triplet_mpi_type,&
    !      & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
    ! CALL MPI_File_read_all(mpi_file_handler, triplet_list%data, local_triplets,&
    !      & triplet_mpi_type, mpi_status,ierr)
    ! CALL MPI_File_close(mpi_file_handler,ierr)
    !
    ! CALL this%FillFromTripletList(triplet_list)
    !
    ! !! Cleanup
    ! CALL triplet_list%Destruct
  END SUBROUTINE ConstructFromBinary_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a distributed sparse matrix to a matrix market file.
  !! Read \cite boisvert1996matrix for the details.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteMatrixToMatrixMarket_psr(this,file_name)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
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
    CALL this%MergeBlocks(merged_local_data)

    bytes_per_character = sizeof(temp_char)

    !! Create the matrix size line
    NEW_LINE_LENGTH = LEN(new_line('A'))
    WRITE(temp_string1,'(A)') "%%MatrixMarket matrix coordinate real general" &
         & //new_line('A')//"%"//new_line('A')
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string1)) :: header_line1)
    header_line1 = TRIM(temp_string1)

    WRITE(temp_string2,*) this%actual_matrix_dimension, &
         & this%actual_matrix_dimension, this%GetSize()
    !! I don't understand why the +1 is needed, but it is.
    ALLOCATE(CHARACTER(&
         & len=LEN_TRIM(temp_string2)+NEW_LINE_LENGTH+1) :: header_line2)
    WRITE(header_line2,*) TRIM(temp_string2)//new_line('A')

    header_size = LEN(header_line1) + LEN(header_line2)

    !! Local Data
    CALL merged_local_data%ConvertToTripletList(triplet_list)

    !! Absolute Positions
    CALL triplet_list%Shift(this%start_row - 1, this%start_column - 1)

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
  END SUBROUTINE WriteMatrixToMatrixMarket_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a binary file.
  !! Faster than text, so this is good for check pointing.
  !! @param[in] this the Matrix to write.
  !! @param[in] file_name name of the file to write to.
  SUBROUTINE WriteMatrixToBinary_psr(this,file_name)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    CHARACTER(len=*), INTENT(IN) :: file_name
    ! !! Local Data
    ! TYPE(TripletList_r) :: triplet_list
    ! INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
    ! INTEGER :: mpi_file_handler
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: local_data_size
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: write_offset
    ! !! For The Special Datatype
    ! INTEGER :: triplet_mpi_type
    ! !! Temporary Variables
    ! INTEGER :: temp_int
    ! REAL(NTREAL) :: temp_double
    ! INTEGER :: bytes_per_int, bytes_per_double
    ! INTEGER, DIMENSION(3) :: header_buffer
    ! INTEGER :: mpi_status(MPI_STATUS_SIZE)
    ! INTEGER(KIND=MPI_OFFSET_KIND) :: zero_offset = 0
    ! INTEGER :: counter
    ! TYPE(Matrix_lsr) :: merged_local_data
    ! INTEGER :: ierr
    ! TYPE(Triplet_r) :: temp_triplet
    !
    ! !! Merge all the local data
    ! CALL this%MergeBlocks(merged_local_data)
    !
    ! !! Determine Write Location
    ! bytes_per_int = sizeof(temp_int)
    ! bytes_per_double = sizeof(temp_double)
    ! local_data_size = SIZE(merged_local_data%values)*(bytes_per_int*2 + &
    !      & bytes_per_double*1)
    ! header_size = bytes_per_int*3
    ! ALLOCATE(local_values_buffer(this%process_grid%slice_size))
    ! CALL MPI_Allgather(SIZE(merged_local_data%values),1,MPI_INT,&
    !      & local_values_buffer,1,MPI_INT,&
    !      & this%process_grid%within_slice_comm,ierr)
    ! write_offset = 0
    ! write_offset = write_offset + header_size
    ! DO counter = 1,this%process_grid%within_slice_rank
    !    write_offset = write_offset + &
    !         & local_values_buffer(counter)*(bytes_per_int*2+bytes_per_double*1)
    ! END DO
    !
    ! !! Write The File
    ! IF (this%process_grid%between_slice_rank .EQ. 0) THEN
    !    !! Create Special MPI Type
    !    triplet_mpi_type = temp_triplet%GetMPIType()
    !
    !    CALL merged_local_data%ConvertToTripletList(triplet_list)
    !    !! Absolute Positions
    !    CALL triplet_list%Shift(this%start_row -1, this%start_column-1)
    !    CALL MPI_File_open(this%process_grid%within_slice_comm,file_name,&
    !         & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY),MPI_INFO_NULL, &
    !         & mpi_file_handler, ierr)
    !    !! Write Header
    !    IF (this%process_grid%within_slice_rank .EQ. 0) THEN
    !       header_buffer(1) = this%actual_matrix_dimension
    !       header_buffer(2) = this%actual_matrix_dimension
    !       header_buffer(3) = SUM(local_values_buffer)
    !       CALL MPI_File_write_at(mpi_file_handler,zero_offset,header_buffer,3,&
    !            & MPI_INT,mpi_status,ierr)
    !    END IF
    !    !! Write The Rest
    !    CALL MPI_File_set_view(mpi_file_handler,write_offset,triplet_mpi_type,&
    !         & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
    !    CALL MPI_File_write(mpi_file_handler,triplet_list%data, &
    !         & triplet_list%CurrentSize, triplet_mpi_type,MPI_STATUS_IGNORE, &
    !         & ierr)
    !
    !    !! Cleanup
    !    CALL MPI_File_close(mpi_file_handler,ierr)
    !    CALL triplet_list%Destruct
    ! END IF
    ! DEALLOCATE(local_values_buffer)
    ! CALL MPI_Barrier(this%process_grid%global_comm,ierr)
    ! CALL merged_local_data%Destruct
  END SUBROUTINE WriteMatrixToBinary_psr
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
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CLASS(TripletList), INTENT(IN) :: triplet_list
    LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    !! Local Data
    TYPE(Permutation_t) :: basic_permutation
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix
    TYPE(Matrix_lsr) :: gathered_matrix
    TYPE(ReduceHelper_t) :: gather_helper
    REAL(NTREAL), PARAMETER :: threshold = 0.0
    LOGICAL :: preduplicated

    SELECT TYPE(triplet_list)
    CLASS IS (TripletList_r)

    IF (.NOT. PRESENT(preduplicated_in)) THEN
       preduplicated = .FALSE.
    ELSE
       preduplicated = preduplicated_in
    END IF

    !! First we redistribute the triplet list to get all the local data
    !! on the correct process.
    CALL basic_permutation%InitDefault(this%logical_matrix_dimension)
    CALL RedistributeData(this,basic_permutation%index_lookup, &
         & basic_permutation%reverse_index_lookup, triplet_list, &
         & sorted_triplet_list)

    !! Now we can just construct a local matrix.
    CALL local_matrix%InitEmpty(this%local_rows, this%local_columns)
    !! And reduce over the Z dimension. This can be accomplished by
    !! summing up.
    IF (.NOT. PRESENT(preduplicated_in) .OR. .NOT. preduplicated_in) THEN
       CALL ReduceMatrixSizes(local_matrix, &
            & this%process_grid%between_slice_comm, gather_helper)
       DO WHILE(.NOT. gather_helper%TestSize())
       END DO
       CALL ReduceAndSumMatrixData(local_matrix, gathered_matrix, &
            & this%process_grid%between_slice_comm, gather_helper)
       DO WHILE(.NOT. gather_helper%TestOuter())
       END DO
       DO WHILE(.NOT. gather_helper%TestInner())
       END DO
       DO WHILE(.NOT. gather_helper%TestData())
       END DO
       CALL ReduceAndSumMatrixCleanup(local_matrix, gathered_matrix, threshold,&
            & gather_helper)
       CALL this%SplitToLocal(gathered_matrix)
    ELSE
       CALL this%SplitToLocal(local_matrix)
    END IF

    CALL local_matrix%Destruct
    CALL sorted_triplet_list%Destruct

    END SELECT
  END SUBROUTINE FillMatrixFromTripletList_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  !! @param[inout] this the matrix being filled.
  SUBROUTINE FillMatrixIdentity_psr(this)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    INTEGER :: i, j
    INTEGER :: total_values
    TYPE(Matrix_lsr) :: local_matrix

    !! There can't be more than one entry per row
    CALL triplet_list%Init(this%local_rows)

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
    CALL unsorted_triplet_list%Init(total_values)
    unsorted_triplet_list%data = triplet_list%data(:total_values)
    CALL SortTripletList(unsorted_triplet_list,this%local_columns,&
         & this%local_rows, sorted_triplet_list)
    CALL local_matrix%InitFromTripletList(sorted_triplet_list, &
         & this%local_rows, this%local_columns)

    CALL this%SplitToLocal(local_matrix)

    CALL local_matrix%Destruct
    CALL triplet_list%Destruct
    CALL unsorted_triplet_list%Destruct
    CALL sorted_triplet_list%Destruct
  END SUBROUTINE FillMatrixIdentity_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !! If you don't specify permuterows, will default to permuting rows.
  !! @param[inout] this the matrix being filled.
  !! @param[in] permutation_vector describes for each row/column, where it goes.
  !! @param[in] permuterows if true permute rows, false permute columns.
  SUBROUTINE FillMatrixPermutation_psr(this, permutation_vector, &
       & permuterows)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    LOGICAL, OPTIONAL, INTENT(IN) :: permuterows
    !! Local Data
    LOGICAL :: rows
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    INTEGER :: total_values
    INTEGER :: counter
    INTEGER :: local_row, local_column
    TYPE(Matrix_lsr) :: local_matrix

    !! Figure out what type of permutation
    IF (PRESENT(permuterows) .AND. permuterows .EQV. .FALSE.) THEN
       rows = .FALSE.
    ELSE
       rows = .TRUE.
    END IF

    !! Build Local Triplet List
    !! There can't be more than one entry per row
    CALL triplet_list%Init(this%local_rows)
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
    CALL unsorted_triplet_list%Init(total_values)
    unsorted_triplet_list%data = triplet_list%data(:total_values)
    CALL SortTripletList(unsorted_triplet_list, this%local_columns, &
         & this%local_rows, sorted_triplet_list)
    CALL local_matrix%InitFromTripletList(sorted_triplet_list, &
         & this%local_rows, this%local_columns)

    CALL this%SplitToLocal(local_matrix)

    CALL local_matrix%Destruct
    CALL triplet_list%Destruct
    CALL unsorted_triplet_list%Destruct
    CALL sorted_triplet_list%Destruct
  END SUBROUTINE FillMatrixPermutation_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !! Data is returned with absolute coordinates.
  !! @param[in] this the distributed sparse matrix.
  !! @param[inout] triplet_list the list to fill.
  PURE SUBROUTINE GetTripletList_psr(this, triplet_list)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    CLASS(TripletList), INTENT(INOUT) :: triplet_list
    !! Local Data
    TYPE(Matrix_lsr) :: merged_local_data

    SELECT TYPE(triplet_list)
    CLASS IS (TripletList_r)
    !! Merge all the local data
    CALL this%MergeBlocks(merged_local_data)

    CALL merged_local_data%ConvertToTripletList(triplet_list)
    CALL triplet_list%Shift(this%start_row-1, this%start_column-1)
    END SELECT
  END SUBROUTINE GetTripletList_psr
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
  SUBROUTINE GetBlock_psr(this, triplet_list, start_row, end_row, &
       & start_column, end_column)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    CLASS(TripletList), INTENT(INOUT) :: triplet_list
    INTEGER, INTENT(IN) :: start_row, end_row
    INTEGER, INTENT(IN) :: start_column, end_column
    !! Local Data
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(TripletList_r) :: local_triplet_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_end_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_end_list
    !! Send Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    !! Receive Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! Temporary
    INTEGER :: counter, p_counter
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: ierr

    SELECT TYPE(triplet_list)
    CLASS IS (TripletList_r)

    !! Merge all the local data
    CALL this%MergeBlocks(merged_local_data)
    CALL merged_local_data%ConvertToTripletList(local_triplet_list)

    !! Share the start row/column information across processes
    ALLOCATE(row_start_list(this%process_grid%slice_size))
    ALLOCATE(column_start_list(this%process_grid%slice_size))
    ALLOCATE(row_end_list(this%process_grid%slice_size))
    ALLOCATE(column_end_list(this%process_grid%slice_size))
    CALL MPI_Allgather(start_row, 1, MPI_INT, row_start_list, 1, MPI_INT, &
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allgather(start_column, 1, MPI_INT, column_start_list, 1, MPI_INT,&
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allgather(end_row, 1, MPI_INT, row_end_list, 1, MPI_INT, &
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allgather(end_column, 1, MPI_INT, column_end_list, 1, MPI_INT,&
         & this%process_grid%within_slice_comm, ierr)

    !! Count The Number of Elements To Send To Each Process
    ALLOCATE(send_per_proc(this%process_grid%slice_size))
    send_per_proc = 0
    DO counter = 1, local_triplet_list%CurrentSize
       CALL local_triplet_list%Get(counter, temp_triplet)
       temp_triplet%index_row = temp_triplet%index_row + this%start_row - 1
       temp_triplet%index_column = temp_triplet%index_column + &
            & this%start_column - 1
       DO p_counter = 1, this%process_grid%slice_size
          IF (temp_triplet%index_row .GE. row_start_list(p_counter) .AND. &
               & temp_triplet%index_row .LT. row_end_list(p_counter) .AND. &
               & temp_triplet%index_column .GE. column_start_list(p_counter) .AND. &
               & temp_triplet%index_column .LT. column_end_list(p_counter)) THEN
             send_per_proc(p_counter) = send_per_proc(p_counter) + 1
             EXIT
          END IF
       END DO
    END DO
    !! Compute send buffer offsets
    ALLOCATE(send_buffer_offsets(this%process_grid%slice_size))
    send_buffer_offsets(1) = 1
    DO counter = 2, this%process_grid%slice_size
       send_buffer_offsets(counter) = send_buffer_offsets(counter-1) + &
            & send_per_proc(counter-1)
    END DO

    !! Build a send buffer
    ALLOCATE(send_buffer_row(local_triplet_list%CurrentSize))
    ALLOCATE(send_buffer_col(local_triplet_list%CurrentSize))
    ALLOCATE(send_buffer_val(local_triplet_list%CurrentSize))
    DO counter = 1, local_triplet_list%CurrentSize
       CALL local_triplet_list%Get(counter, temp_triplet)
       temp_triplet%index_row = temp_triplet%index_row + this%start_row - 1
       temp_triplet%index_column = temp_triplet%index_column + &
            & this%start_column - 1
       DO p_counter = 1, this%process_grid%slice_size
          IF (temp_triplet%index_row .GE. row_start_list(p_counter) .AND. &
               & temp_triplet%index_row .LT. row_end_list(p_counter) .AND. &
               & temp_triplet%index_column .GE. column_start_list(p_counter) .AND. &
               & temp_triplet%index_column .LT. column_end_list(p_counter)) THEN
             send_buffer_row(send_buffer_offsets(p_counter)) = &
                  & temp_triplet%index_row
             send_buffer_col(send_buffer_offsets(p_counter)) = &
                  & temp_triplet%index_column
             send_buffer_val(send_buffer_offsets(p_counter)) = &
                  & temp_triplet%point_value
             send_buffer_offsets(p_counter) = send_buffer_offsets(p_counter) + 1
             EXIT
          END IF
       END DO
    END DO
    !! Reset send buffer offsets. But since we're using MPI now, use zero
    !! based indexing.
    send_buffer_offsets(1) = 0
    DO counter = 2, this%process_grid%slice_size
       send_buffer_offsets(counter) = send_buffer_offsets(counter-1) + &
            & send_per_proc(counter-1)
    END DO

    !! Build a receive buffer
    ALLOCATE(recv_per_proc(this%process_grid%slice_size))
    CALL MPI_Alltoall(send_per_proc, 1, MPI_INT, recv_per_proc, 1, MPI_INT, &
         & this%process_grid%within_slice_comm, ierr)
    ALLOCATE(recv_buffer_offsets(this%process_grid%slice_size))
    recv_buffer_offsets(1) = 0
    DO counter = 2, this%process_grid%slice_size
       recv_buffer_offsets(counter) = recv_buffer_offsets(counter-1) + &
            & recv_per_proc(counter-1)
    END DO
    ALLOCATE(recv_buffer_row(SUM(recv_per_proc)))
    ALLOCATE(recv_buffer_col(SUM(recv_per_proc)))
    ALLOCATE(recv_buffer_val(SUM(recv_per_proc)))

    !! Send
    CALL MPI_Alltoallv(send_buffer_row, send_per_proc, send_buffer_offsets, &
         & MPI_INT, recv_buffer_row, recv_per_proc, recv_buffer_offsets, &
         & MPI_INT, this%process_grid%within_slice_comm, ierr)
    CALL MPI_Alltoallv(send_buffer_col, send_per_proc, send_buffer_offsets, &
         & MPI_INT, recv_buffer_col, recv_per_proc, recv_buffer_offsets, &
         & MPI_INT, this%process_grid%within_slice_comm, ierr)
    CALL MPI_Alltoallv(send_buffer_val, send_per_proc, send_buffer_offsets, &
         & MPINTREAL, recv_buffer_val, recv_per_proc, recv_buffer_offsets, &
         & MPINTREAL, this%process_grid%within_slice_comm, ierr)

    !! Convert receive buffer to triplet list
    CALL triplet_list%Init(SUM(recv_per_proc))
    DO counter=1, SUM(recv_per_proc)
       triplet_list%data(counter)%index_row = recv_buffer_row(counter)
       triplet_list%data(counter)%index_column = recv_buffer_col(counter)
       triplet_list%data(counter)%point_value = recv_buffer_val(counter)
    END DO

    !! Cleanup
    IF (ALLOCATED(row_start_list)) DEALLOCATE(row_start_list)
    IF (ALLOCATED(column_start_list)) DEALLOCATE(column_start_list)
    IF (ALLOCATED(row_end_list)) DEALLOCATE(row_end_list)
    IF (ALLOCATED(column_end_list)) DEALLOCATE(column_end_list)
    IF (ALLOCATED(recv_buffer_offsets)) DEALLOCATE(recv_buffer_offsets)
    IF (ALLOCATED(recv_buffer_val)) DEALLOCATE(recv_buffer_val)
    IF (ALLOCATED(recv_buffer_col)) DEALLOCATE(recv_buffer_col)
    IF (ALLOCATED(recv_buffer_row)) DEALLOCATE(recv_buffer_row)
    IF (ALLOCATED(recv_per_proc)) DEALLOCATE(recv_per_proc)
    IF (ALLOCATED(send_buffer_val)) DEALLOCATE(send_buffer_val)
    IF (ALLOCATED(send_buffer_col)) DEALLOCATE(send_buffer_col)
    IF (ALLOCATED(send_buffer_row)) DEALLOCATE(send_buffer_row)
    IF (ALLOCATED(send_buffer_offsets)) DEALLOCATE(send_buffer_offsets)
    IF (ALLOCATED(send_per_proc)) DEALLOCATE(send_per_proc)

    END SELECT
  END SUBROUTINE GetBlock_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out information about a distributed sparse matrix.
  !! Sparsity, and load balancing information.
  !! @param[in] this the matrix to print information about.
  SUBROUTINE PrintMatrixInformation_psr(this)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    CALL this%GetLoadBalance(min_size, max_size)
    sparsity = REAL(this%GetSize(),KIND=NTREAL) / &
         & (REAL(this%actual_matrix_dimension,KIND=NTREAL)**2)

    CALL WriteHeader("Load_Balance")
    CALL EnterSubLog
    CALL WriteListElement(key="min_size", int_value_in=min_size)
    CALL WriteListElement(key="max_size", int_value_in=max_size)
    CALL ExitSubLog
    CALL WriteElement(key="Dimension",int_value_in=this%actual_matrix_dimension)
    CALL WriteElement(key="Sparsity", float_value_in=sparsity)
  END SUBROUTINE PrintMatrixInformation_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a distributed sparse matrix.
  !! This is a serial print routine, and should probably only be used for debug
  !! purposes.
  !! @param[in] this the matrix to print.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintMatrix_psr(this, file_name_in)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Helpers For Communication
    TYPE(ReduceHelper_t) :: row_helper
    TYPE(ReduceHelper_t) :: column_helper
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    !! Temporary Variables
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(Matrix_lsr) :: merged_local_dataT
    TYPE(Matrix_lsr) :: merged_columns
    TYPE(Matrix_lsr) :: merged_columnsT
    TYPE(Matrix_lsr) :: full_gathered
    INTEGER :: ierr

    !! Merge all the local data
    CALL this%MergeBlocks(merged_local_data)

    !! Merge Columns
    CALL merged_local_dataT%Transpose(merged_local_data)
    CALL ReduceMatrixSizes(merged_local_dataT, this%process_grid%column_comm, &
         & column_helper)
    CALL MPI_Wait(column_helper%size_request ,mpi_status,ierr)
    CALL ReduceAndComposeMatrixData(merged_local_dataT, &
         & this%process_grid%column_comm,merged_columns, &
         & column_helper)
    CALL MPI_Wait(column_helper%outer_request,mpi_status,ierr)
    CALL MPI_Wait(column_helper%inner_request,mpi_status,ierr)
    CALL MPI_Wait(column_helper%data_request,mpi_status,ierr)
    CALL ReduceAndComposeMatrixCleanup(merged_local_dataT,merged_columns, &
         & column_helper)

    !! Merge Rows
    CALL merged_columnsT%Transpose(merged_columns)
    CALL ReduceMatrixSizes(merged_columnsT, this%process_grid%row_comm, &
         & row_helper)
    CALL MPI_Wait(row_helper%size_request,mpi_status,ierr)
    CALL ReduceAndComposeMatrixData(merged_columnsT, &
         & this%process_grid%row_comm, full_gathered, row_helper)
    CALL MPI_Wait(row_helper%outer_request,mpi_status,ierr)
    CALL MPI_Wait(row_helper%inner_request,mpi_status,ierr)
    CALL MPI_Wait(row_helper%data_request,mpi_status,ierr)
    CALL ReduceAndComposeMatrixCleanup(merged_columnsT,full_gathered,row_helper)

    !! Make these changes so that it prints the logical rows/columns
    full_gathered%rows = this%actual_matrix_dimension
    full_gathered%columns = this%actual_matrix_dimension

    IF (this%process_grid%IsRoot()) THEN
       IF (PRESENT(file_name_in)) THEN
          CALL full_gathered%Print(file_name_in)
       ELSE
          CALL full_gathered%Print
       END IF
    END IF

    CALL merged_local_data%Destruct
    CALL merged_local_dataT%Destruct
    CALL merged_columns%Destruct
    CALL merged_columnsT%Destruct
  END SUBROUTINE PrintMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get a measure of how load balanced this matrix is. For each process, the
  !! number of non-zero entries is calculated. Then, this function returns
  !! the max and min of those values.
  !! @param[in] this The matrix to compute the measure on.
  !! @param[out] min_size the minimum entries contained on a single process.
  !! @param[out] max_size the maximum entries contained on a single process.
  SUBROUTINE GetMatrixLoadBalance_psr(this, min_size, max_size)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: min_size
    INTEGER, INTENT(OUT) :: max_size
    !! Local Data
    INTEGER :: local_size
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL this%MergeBlocks(merged_local_data)

    local_size = SIZE(merged_local_data%values)
    CALL MPI_Allreduce(local_size,max_size,1,MPI_INT,MPI_MAX,&
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allreduce(local_size,min_size,1,MPI_INT,MPI_MIN,&
         & this%process_grid%within_slice_comm, ierr)

    CALL merged_local_data%Destruct
  END SUBROUTINE GetMatrixLoadBalance_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the total number of non zero entries in the distributed sparse matrix.
  !! @param[in] this the distributed sparse matrix to calculate the non-zero
  !! entries of.
  !! @return the number of non-zero entries in the matrix.
  FUNCTION GetMatrixSize_psr(this) RESULT(total_size)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    INTEGER(c_long) :: total_size
    !! Local Data
    !integer :: local_size
    REAL(NTREAL) :: local_size
    REAL(NTREAL) :: temp_size
    TYPE(Matrix_lsr) :: merged_local_data
    INTEGER :: ierr

    !! Merge all the local data
    CALL this%MergeBlocks(merged_local_data)

    local_size = SIZE(merged_local_data%values)
    CALL MPI_Allreduce(local_size,temp_size,1,MPINTREAL,MPI_SUM,&
         & this%process_grid%within_slice_comm, ierr)

    total_size = INT(temp_size,kind=c_long)

    CALL merged_local_data%Destruct
  END FUNCTION GetMatrixSize_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !! All (absolute) values below the threshold are set to zero.
  !! @param[inout] this matrix to filter
  !! @param[in] threshold (absolute) values below this are filtered
  SUBROUTINE FilterMatrix_psr(this, threshold)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary
    INTEGER :: counter
    INTEGER :: size_temp
    TYPE(ProcessGrid_t) :: grid_temp

    CALL this%GetTripletList(triplet_list)
    CALL new_list%Init
    DO counter=1,triplet_list%CurrentSize
       CALL triplet_list%Get(counter, temporary)
       IF (ABS(temporary%point_value) .GT. threshold) THEN
          CALL new_list%Append(temporary)
       END IF
    END DO
    size_temp = this%actual_matrix_dimension
    grid_temp = this%process_grid
    CALL this%Destruct
    CALL this%InitEmpty(size_temp, grid_temp)
    CALL this%FillFromTripletList(new_list)
  END SUBROUTINE FilterMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge together the local matrix blocks into one big matrix.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[inout] merged_matrix the merged matrix.
  PURE SUBROUTINE MergeMatrixLocalBlocks_psr(this, merged_matrix)
    !! Parameters
    CLASS(Matrix_psr), INTENT(IN) :: this
    CLASS(Matrix_ls), INTENT(INOUT) :: merged_matrix

    SELECT TYPE(merged_matrix)
    CLASS IS (Matrix_lsr)
    CALL ComposeMatrix(this%local_data, &
         & this%process_grid%number_of_blocks_rows, &
         & this%process_grid%number_of_blocks_columns, merged_matrix)
    END SELECT
  END SUBROUTINE MergeMatrixLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[in] matrix_to_split the matrix to split up.
  SUBROUTINE SplitMatrixToLocalBlocks_psr(this, matrix_to_split)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CLASS(Matrix_ls), INTENT(IN) :: matrix_to_split

    SELECT TYPE(matrix_to_split)
    CLASS IS (Matrix_lsr)
    CALL SplitMatrix(matrix_to_split, &
         & this%process_grid%number_of_blocks_rows, &
         & this%process_grid%number_of_blocks_columns, this%local_data)
    END SELECT
  END SUBROUTINE SplitMatrixToLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix.
  !! @param[in] matA The matrix to transpose.
  !! @param[inout] this = matA^T
  SUBROUTINE TransposeMatrix_psr(this, matA)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CLASS(Matrix_p), INTENT(IN) :: matA
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary, temporary_t
    INTEGER :: counter

    CALL new_list%Init
    CALL matA%GetTripletList(triplet_list)
    DO counter=1,triplet_list%CurrentSize
       IF (MOD(counter, matA%process_grid%num_process_slices) .EQ. &
            & matA%process_grid%my_slice) THEN
          CALL triplet_list%Get(counter, temporary)
          temporary_t%index_row = temporary%index_column
          temporary_t%index_column = temporary%index_row
          temporary_t%point_value = temporary%point_value
          CALL new_list%Append(temporary_t)
       END IF
    END DO

    CALL this%Destruct
    CALL this%InitEmpty(matA%actual_matrix_dimension, matA%process_grid)
    CALL this%FillFromTripletList(new_list)
    CALL new_list%Destruct
    CALL triplet_list%Destruct
  END SUBROUTINE TransposeMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split the current communicator, and give each group a complete copy of this
  !! @param[in] this the matrix to split.
  !! @param[out] split_mat a copy of the matrix hosted on a small process grid.
  !! @param[out] my_color distinguishes between the two groups (optional).
  !! @param[out] split_slice if we split along the slice direction, this is True
  SUBROUTINE CommSplitMatrix_psr(this, split_mat, my_color, split_slice)
    !! Parameters
    CLASS(Matrix_psr), INTENT(INOUT) :: this
    CLASS(Matrix_p), INTENT(INOUT) :: split_mat
    INTEGER, INTENT(OUT) :: my_color
    LOGICAL, INTENT(OUT) :: split_slice
    !! For Grid Splitting
    TYPE(ProcessGrid_t) :: new_grid
    INTEGER :: between_grid_comm
    INTEGER :: between_grid_size
    INTEGER :: between_grid_rank
    !! For Data Redistribution
    TYPE(TripletList_r) :: full_list, new_list
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list
    INTEGER :: fsize
    INTEGER :: counter
    INTEGER :: ierr

    IF (this%process_grid%total_processors .EQ. 1) THEN
       CALL split_mat%Copy(this)
       my_color = 0
       split_slice = .TRUE.
    ELSE
       !! Split The Grid
       call this%process_grid%Split(new_grid, my_color, split_slice, &
            & between_grid_comm)

       !! Copy The Data Across New Process Grids. Unnecessary if we just split
       !! by slices.
       CALL this%GetTripletList(full_list)
       IF (.NOT. split_slice) THEN
          CALL MPI_COMM_SIZE(between_grid_comm, between_grid_size, ierr)
          CALL MPI_COMM_RANK(between_grid_comm, between_grid_rank, ierr)

          !! Build Send Lists
          fsize = full_list%CurrentSize
          ALLOCATE(send_list(between_grid_size))
          IF (my_color .EQ. 0) THEN
             !! The smaller process grid only needs to send to process 2
             CALL send_list(1)%Init()
             CALL send_list(2)%init(full_list%CurrentSize)
             send_list(2)%data(:fsize) = full_list%data(:fsize)
             DO counter = 3, between_grid_size
                CALL send_list(counter)%Init()
             END DO
          ELSE
             !! The larger process grid only needs to send to process 1
             CALL send_list(1)%Init(full_list%CurrentSize)
             send_list(1)%data(:fsize) = full_list%data(:fsize)
             DO counter = 2, between_grid_size
                CALL send_list(counter)%Init
             END DO
          END IF
          CALL send_list(between_grid_rank+1)%Init(full_list%CurrentSize)
          send_list(between_grid_rank+1)%data(:fsize) = full_list%data(:fsize)
          CALL RedistributeTripletLists(send_list, between_grid_comm, new_list)
       END IF

       !! Create The New Matrix
       CALL split_mat%InitEmpty(this%actual_matrix_dimension, new_grid)
       IF (.NOT. split_slice) THEN
          CALL split_mat%FillFromTripletList(new_list, .TRUE.)
       ELSE
          CALL split_mat%FillFromTripletList(full_list, .TRUE.)
       END IF

       !! Cleanup
       CALL full_list%Destruct
       CALL new_list%Destruct
       IF (ALLOCATED(send_list)) THEN
          DO counter = 1, between_grid_size
             CALL send_list(counter)%Destruct
          END DO
       END IF
    END IF

  END SUBROUTINE CommSplitMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate a matrix size that can be divided by the number of processors.
  !! @param[in] matrix_dim the dimension of the actual matrix.
  !! @return a new dimension which includes padding.
  !! @todo write a more optimal algorithm using either plain brute force,
  !! or a prime factor based algorithm.
  PURE FUNCTION CalculateScaledDimension(process_grid, matrix_dim) &
       & RESULT(scaled_dim)
    !! Parameters
    TYPE(ProcessGrid_t), INTENT(IN) :: process_grid
    INTEGER, INTENT(IN) :: matrix_dim
    INTEGER :: scaled_dim
    !! Local Data
    INTEGER :: size_ratio
    INTEGER :: lcm

    lcm = process_grid%block_multiplier*process_grid%num_process_slices* &
        & process_grid%num_process_columns*process_grid%num_process_rows

    size_ratio = matrix_dim/lcm
    IF (size_ratio * lcm .EQ. matrix_dim) THEN
       scaled_dim = matrix_dim
    ELSE
       scaled_dim = (size_ratio + 1)*(lcm)
    END IF
  END FUNCTION CalculateScaledDimension
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
  SUBROUTINE RedistributeData(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !! Parameters
    TYPE(Matrix_psr), INTENT(INOUT) :: this
    INTEGER, DIMENSION(:), INTENT(IN) :: index_lookup
    INTEGER, DIMENSION(:), INTENT(IN) :: reverse_index_lookup
    TYPE(TripletList_r), INTENT(IN) :: initial_triplet_list
    TYPE(TripletList_r), INTENT(OUT) :: sorted_triplet_list
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_lookup
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_lookup
    INTEGER, DIMENSION(:), ALLOCATABLE :: location_list_within_slice
    TYPE(TripletList_r) :: gathered_list
    TYPE(TripletList_r), DIMENSION(this%process_grid%slice_size) :: send_triplet_lists
    !! Temporary Values
    INTEGER :: row_size, column_size
    INTEGER :: temp_row, temp_column
    INTEGER :: process_id
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: counter

    !! First we need to figure out where our local elements go
    ALLOCATE(row_lookup(SIZE(index_lookup)))
    ALLOCATE(column_lookup(SIZE(index_lookup)))
    row_size = SIZE(index_lookup)/this%process_grid%num_process_rows
    DO counter = LBOUND(index_lookup,1), UBOUND(index_lookup,1)
       row_lookup(index_lookup(counter)) = (counter-1)/(row_size)
    END DO
    column_size = SIZE(index_lookup)/this%process_grid%num_process_columns
    DO counter = LBOUND(index_lookup,1), UBOUND(index_lookup,1)
       column_lookup(index_lookup(counter)) = (counter-1)/(column_size)
    END DO
    ALLOCATE(location_list_within_slice(initial_triplet_list%CurrentSize))
    DO counter = 1, initial_triplet_list%CurrentSize
       temp_row = row_lookup(initial_triplet_list%data(counter)%index_row)
       temp_column = &
            & column_lookup(initial_triplet_list%data(counter)%index_column)
       location_list_within_slice(counter) = &
            & temp_column+temp_row*this%process_grid%num_process_columns
    END DO

    !! Build A Send Buffer
    DO counter = 1, this%process_grid%slice_size
       CALL send_triplet_lists(counter)%Init
    END DO
    DO counter = 1, initial_triplet_list%CurrentSize
       process_id = location_list_within_slice(counter)
       CALL initial_triplet_list%Get(counter, temp_triplet)
       CALL send_triplet_lists(process_id+1)%Append(temp_triplet)
    END DO

    !! Actual Send
    CALL RedistributeTripletLists(send_triplet_lists, &
         & this%process_grid%within_slice_comm, gathered_list)

    !! Adjust Indices to Local
    DO counter = 1, gathered_list%CurrentSize
       gathered_list%data(counter)%index_row = &
            & reverse_index_lookup(gathered_list%data(counter)%index_row) - &
            & this%start_row + 1
       gathered_list%data(counter)%index_column = &
            & reverse_index_lookup(gathered_list%data(counter)%index_column) - &
            & this%start_column + 1
    END DO
    CALL SortTripletList(gathered_list, this%local_columns, this%local_rows, &
         & sorted_triplet_list)

    !! Cleanup
    DO counter = 1, this%process_grid%slice_size
       CALL send_triplet_lists(counter)%Destruct
    END DO
    DEALLOCATE(row_lookup)
    DEALLOCATE(column_lookup)
    DEALLOCATE(location_list_within_slice)
    CALL gathered_list%Destruct
  END SUBROUTINE RedistributeData
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixPSModule
