!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Operations.
MODULE PSMatrixModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX, &
       & MPINTINTEGER
  USE ErrorModule, ONLY : Error_t, ConstructError, SetGenericError, &
       & CheckMPIError
  USE LoggingModule, ONLY : &
       & EnterSubLog, ExitSubLog, WriteElement, WriteListElement, WriteHeader
  USE MatrixMarketModule, ONLY : ParseMMHeader, MM_COMPLEX
  USE PermutationModule, ONLY : Permutation_t, ConstructDefaultPermutation
  USE ProcessGridModule, ONLY : ProcessGrid_t, global_grid, IsRoot, &
       & SplitProcessGrid
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & ReduceANdSumMatrixSizes, ReduceAndSumMatrixData, &
       & ReduceAndSumMatrixCleanup, TestReduceSizeRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, DestructMatrix, &
       & PrintMatrix, TransposeMatrix, ConjugateMatrix, SplitMatrix, &
       & ComposeMatrix, ConvertMatrixType, MatrixToTripletList, &
       & ConstructMatrixFromTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletModule, ONLY : Triplet_r, Triplet_c, GetMPITripletType_r, &
       & GetMPITripletType_c
  USE TripletListModule, ONLY : TripletList_r, TripletList_c, &
       & ConstructTripletList, &
       & DestructTripletList, SortTripletList, AppendToTripletList, &
       & SymmetrizeTripletList, GetTripletAt, RedistributeTripletLists, &
       & ShiftTripletList
  USE NTMPIModule
  USE, INTRINSIC :: ISO_C_BINDING
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
     TYPE(ProcessGrid_t), POINTER :: process_grid !< process grid to operate on
     LOGICAL :: is_complex !< true if the matrix data is true.
  END TYPE Matrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Constructors/Destructors
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: DestructMatrix
  PUBLIC :: CopyMatrix
  PUBLIC :: SetMatrixProcessGrid
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
  PUBLIC :: GetMatrixSlice
  !! Printing To The Console
  PUBLIC :: PrintMatrix
  PUBLIC :: PrintMatrixInformation
  !! Utilities
  PUBLIC :: ConvertMatrixToReal
  PUBLIC :: ConvertMatrixToComplex
  PUBLIC :: GetMatrixLoadBalance
  PUBLIC :: GetMatrixSize
  PUBLIC :: FilterMatrix
  PUBLIC :: MergeMatrixLocalBlocks
  PUBLIC :: SplitMatrixToLocalBlocks
  PUBLIC :: TransposeMatrix
  PUBLIC :: ConjugateMatrix
  PUBLIC :: CommSplitMatrix
  PUBLIC :: ResizeMatrix
  PUBLIC :: GatherMatrixToProcess
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
  INTERFACE ConjugateMatrix
     MODULE PROCEDURE ConjugateMatrix_ps
  END INTERFACE
  INTERFACE CommSplitMatrix
     MODULE PROCEDURE CommSplitMatrix_ps
  END INTERFACE
  INTERFACE GatherMatrixToProcess
     MODULE PROCEDURE GatherMatrixToProcess_psr
     MODULE PROCEDURE GatherMatrixToProcess_psc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty sparse, distributed, matrix.
  SUBROUTINE ConstructEmptyMatrix_ps(this, matrix_dim_, process_grid_in, &
       & is_complex_in)
    !> The matrix to be constructed.
    TYPE(Matrix_ps), INTENT(INOUT)            :: this
    !> The dimension of the full matrix.
    INTEGER, INTENT(IN)                       :: matrix_dim_
    !> True if you want to use complex numbers.
    LOGICAL, INTENT(IN), OPTIONAL             :: is_complex_in
    !> A process grid to host the matrix.
    TYPE(ProcessGrid_t), INTENT(IN), TARGET, OPTIONAL :: process_grid_in
    !! Local Variables
    TYPE(Matrix_lsr) :: zeromatrix_r
    TYPE(Matrix_lsc) :: zeromatrix_c

    CALL DestructMatrix(this)

    !! Process Grid
    IF (PRESENT(process_grid_in)) THEN
       this%process_grid => process_grid_in
    ELSE
       this%process_grid => global_grid
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
  !> to determine the parameters. Note that no data is copied, the matrix
  !> will be empty.
  SUBROUTINE ConstructEmptyMatrix_ps_cp(this, reference_matrix)
    !! Parameters
    !> The matrix to be constructed.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The reference matrix to take parameters from.
    TYPE(Matrix_ps), INTENT(IN) :: reference_matrix

    CALL ConstructEmptyMatrix(this, reference_matrix%actual_matrix_dimension, &
         & reference_matrix%process_grid, reference_matrix%is_complex)
  END SUBROUTINE ConstructEmptyMatrix_ps_cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a distributed sparse matrix.
  PURE SUBROUTINE DestructMatrix_ps(this)
    !> The matrix to destruct.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: II, JJ

    IF (ALLOCATED(this%local_data_r)) THEN
       DO II = 1, SIZE(this%local_data_r,DIM=1)
          DO JJ = 1, SIZE(this%local_data_r,DIM=2)
             CALL DestructMatrix(this%local_data_r(II,JJ))
          END DO
       END DO
       DEALLOCATE(this%local_data_r)
    END IF

    IF (ALLOCATED(this%local_data_c)) THEN
       DO II = 1, SIZE(this%local_data_c,DIM=1)
          DO JJ = 1, SIZE(this%local_data_c,DIM=2)
             CALL DestructMatrix(this%local_data_c(II,JJ))
          END DO
       END DO
       DEALLOCATE(this%local_data_c)
    END IF

  END SUBROUTINE DestructMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a distributed sparse matrix in a safe way.
  SUBROUTINE CopyMatrix_ps(matA, matB)
    !> The matrix to copy.
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    !> matB = matA.
    TYPE(Matrix_ps), INTENT(INOUT) :: matB

    CALL DestructMatrix(matB)
    matB = matA
  END SUBROUTINE CopyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> When you want to change the process grid of a matrix, you can call
  !> this routine with the new process grid value. Data will be automatically
  !> redistributed.
  SUBROUTINE SetMatrixProcessGrid(this, grid)
    !> The matrix to set the grid of.
    TYPE(Matrix_ps), INTENT(INOUT)  :: this
    !> The grid to set it to.
    TYPE(ProcessGrid_t), INTENT(IN) :: grid
    !! Local variables
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    TYPE(Matrix_ps) :: new_mat

    !! Get the data in a triplet list
    CALL ConstructTripletList(triplet_list_c)
    CALL ConstructTripletList(triplet_list_r)

    IF (this%process_grid%my_slice .EQ. 0) THEN
       IF (this%is_complex) THEN
          CALL GetMatrixTripletList(this, triplet_list_c)
       ELSE
          CALL GetMatrixTripletList(this, triplet_list_r)
       END IF
    END IF

    !! Fill The New Matrix
    CALL ConstructEmptyMatrix(new_mat, this%actual_matrix_dimension, grid, &
         & this%is_complex)
    IF (this%is_complex) THEN
       CALL FillMatrixFromTripletList(new_mat, triplet_list_c)
    ELSE
       CALL FillMatrixFromTripletList(new_mat, triplet_list_r)
    END IF

    !! Copy back to finish
    CALL CopyMatrix(new_mat, this)

    !! Cleanup
    CALL DestructMatrix(new_mat)
    CALL DestructTripletList(triplet_list_c)
    CALL DestructTripletList(triplet_list_r)
  END SUBROUTINE SetMatrixProcessGrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct distributed sparse matrix from a matrix market file in parallel.
  !> Read \cite boisvert1996matrix for the details.
  RECURSIVE SUBROUTINE ConstructMatrixFromMatrixMarket_ps(this, file_name, &
       & process_grid_in)
    !> The file being constructed.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Grid to distribute the matrix on.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: process_grid_in
    !> The name of the file to read.
    CHARACTER(len=*), INTENT(IN) :: file_name
    INTEGER, PARAMETER :: MAX_LINE_LENGTH = 100
    !! File Handles
    INTEGER :: local_file_handler
    INTEGER :: mpi_file_handler
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Reading The File
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(Triplet_r) :: temp_triplet_r
    TYPE(TripletList_c) :: triplet_list_c
    TYPE(Triplet_c) :: temp_triplet_c
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
    REAL(NTREAL) :: realval, cval
    INTEGER :: bytes_per_character
    LOGICAL :: found_comment_line
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER :: full_buffer_counter
    LOGICAL :: end_of_buffer
    LOGICAL :: header_success
    INTEGER :: ierr
    TYPE(Error_t) :: err

    IF (.NOT. PRESENT(process_grid_in)) THEN
       CALL ConstructMatrixFromMatrixMarket(this, file_name, global_grid)
    ELSE
       CALL ConstructError(err)
       !! Setup Involves Just The Root Opening And Reading Parameter Data
       CALL StartTimer("MPI Read Text")
       CALL MPI_Type_size(MPI_CHARACTER, bytes_per_character, ierr)
       IF (IsRoot(process_grid_in)) THEN
          header_length = 0
          local_file_handler = 16
          OPEN(local_file_handler, file=file_name, iostat=ierr, status="old")
          IF (ierr .NE. 0) THEN
             CALL SetGenericError(err, TRIM(file_name)//" doesn't exist", .TRUE.)
          END IF
          !! Parse the header.
          READ(local_file_handler,fmt='(A)') input_buffer
          header_success = ParseMMHeader(input_buffer, sparsity_type, &
               & data_type, pattern_type)
          IF (.NOT. header_success) THEN
             CALL SetGenericError(err, "Invalid File Header", .TRUE.)
          END IF
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
       CALL MPI_Bcast(matrix_rows, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(matrix_columns, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(total_values, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(header_length, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(sparsity_type, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(data_type, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(pattern_type, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)

       !! Build Local Storage
       CALL ConstructEmptyMatrix(this, matrix_rows, process_grid_in, &
            & is_complex_in = (data_type .EQ. MM_COMPLEX))

       !! Global read
       CALL MPI_File_open(this%process_grid%global_comm, file_name, &
            & MPI_MODE_RDONLY, MPI_INFO_NULL,mpi_file_handler,ierr)
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
          IF (local_offset + local_data_size_plus_buffer .GT. &
               & total_file_size) THEN
             local_data_size_plus_buffer = (total_file_size - local_offset)
          END IF
       ELSE
          local_data_size_plus_buffer = 0
       END IF

       !! A buffer to read the data into.
       ALLOCATE(CHARACTER(LEN=local_data_size_plus_buffer) :: mpi_input_buffer)

       !! Do Actual Reading
       CALL MPI_File_read_at_all(mpi_file_handler, local_offset, &
            & mpi_input_buffer, INT(local_data_size_plus_buffer), &
            & MPI_CHARACTER, mpi_status, ierr)

       !! Trim Off The Half Read Line At The Start
       IF (.NOT. this%process_grid%global_rank .EQ. &
            & this%process_grid%RootID) THEN
          full_buffer_counter = INDEX(mpi_input_buffer,new_line('A')) + 1
       ELSE
          full_buffer_counter = 1
       END IF

       !! Read By Line
       end_of_buffer = .FALSE.
       IF (local_data_size_plus_buffer .EQ. 0) THEN
          end_of_buffer = .TRUE.
       END IF

       IF (this%is_complex) THEN
          CALL ConstructTripletList(triplet_list_c)
       ELSE
          CALL ConstructTripletList(triplet_list_r)
       END IF
       DO WHILE(.NOT. end_of_buffer)
          current_line_length = INDEX(mpi_input_buffer(full_buffer_counter:),&
               new_line('A'))

          IF (current_line_length .EQ. 0) THEN !! Hit The End Of The Buffer
             end_of_buffer = .TRUE.
          ELSE
             temp_substring = mpi_input_buffer(full_buffer_counter: &
                  & full_buffer_counter+current_line_length-1)
             IF (current_line_length .GT. 1) THEN
                IF (data_type .EQ. MM_COMPLEX) THEN
                   READ(temp_substring(:current_line_length-1),*) &
                        & temp_triplet_c%index_row, &
                        & temp_triplet_c%index_column, &
                        & realval, cval
                   temp_triplet_c%point_value = &
                        & CMPLX(realval, cval, KIND=NTCOMPLEX)
                   CALL AppendToTripletList(triplet_list_c, temp_triplet_c)
                ELSE
                   READ(temp_substring(:current_line_length-1),*) &
                        & temp_triplet_r%index_row, &
                        & temp_triplet_r%index_column, &
                        & temp_triplet_r%point_value
                   CALL AppendToTripletList(triplet_list_r, temp_triplet_r)
                END IF
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
       IF (this%is_complex) THEN
          CALL SymmetrizeTripletList(triplet_list_c, pattern_type)
          CALL FillMatrixFromTripletList(this,triplet_list_c)
          CALL DestructTripletList(triplet_list_c)
       ELSE
          CALL SymmetrizeTripletList(triplet_list_r, pattern_type)
          CALL FillMatrixFromTripletList(this,triplet_list_r)
          CALL DestructTripletList(triplet_list_r)
       END IF

       DEALLOCATE(mpi_input_buffer)
    END IF

  END SUBROUTINE ConstructMatrixFromMatrixMarket_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a distributed sparse matrix from a binary file in parallel.
  !> Faster than text, so this is good for check pointing.
  RECURSIVE SUBROUTINE ConstructMatrixFromBinary_ps(this, file_name, &
       & process_grid_in)
    !! Parameters
    !> The file being constructed.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Grid to distribute the matrix on.
    TYPE(ProcessGrid_t), INTENT(IN), OPTIONAL :: process_grid_in
    !> The name of the file to read.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! Local Data
    INTEGER :: triplet_mpi_type
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! File Handles
    INTEGER :: mpi_file_handler
    !! Reading The File
    INTEGER :: matrix_rows, matrix_columns, total_values, complex_flag
    INTEGER, DIMENSION(4) :: matrix_information
    INTEGER :: local_triplets
    INTEGER(KIND=MPI_OFFSET_KIND) :: local_offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
    INTEGER :: bytes_per_int, bytes_per_data
    !! Temporary variables
    INTEGER :: mpi_status(MPI_STATUS_SIZE)
    INTEGER :: ierr
    TYPE(Error_t) :: err
    LOGICAL :: error_occured

    IF (.NOT. PRESENT(process_grid_in)) THEN
       CALL ConstructMatrixFromBinary(this, file_name, global_grid)
    ELSE
       CALL ConstructError(err)
       CALL StartTimer("MPI Read Binary")
       CALL MPI_File_open(process_grid_in%global_comm, file_name, &
            & MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_file_handler, ierr)
       error_occured = CheckMPIError(err, TRIM(file_name)//" doesn't exist", &
            & ierr, .TRUE.)

       !! Get The Matrix Parameters
       IF (IsRoot(process_grid_in)) THEN
          local_offset = 0
          CALL MPI_File_read_at(mpi_file_handler, local_offset, &
               & matrix_information, 4, MPINTINTEGER, mpi_status, ierr)
          matrix_rows = matrix_information(1)
          matrix_columns = matrix_information(2)
          total_values = matrix_information(3)
          complex_flag = matrix_information(4)
       END IF

       !! Broadcast Parameters
       CALL MPI_Bcast(matrix_rows, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(matrix_columns, 1, MPINTINTEGER, process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(total_values, 1, MPINTINTEGER ,process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)
       CALL MPI_Bcast(complex_flag, 1, MPINTINTEGER ,process_grid_in%RootID, &
            & process_grid_in%global_comm, ierr)

       !! Build Local Storage
       IF (complex_flag .EQ. 1) THEN
          CALL ConstructEmptyMatrix(this, matrix_rows, process_grid_in, &
               & is_complex_in=.TRUE.)
       ELSE
          CALL ConstructEmptyMatrix(this, matrix_rows, process_grid_in, &
               & is_complex_in=.FALSE.)
       END IF

       CALL MPI_Type_extent(MPINTINTEGER,bytes_per_int,ierr)
       IF (this%is_complex) THEN
          CALL MPI_Type_extent(MPINTCOMPLEX,bytes_per_data,ierr)
          triplet_mpi_type = GetMPITripletType_c()
       ELSE
          CALL MPI_Type_extent(MPINTREAL,bytes_per_data,ierr)
          triplet_mpi_type = GetMPITripletType_r()
       END IF

       !! Compute Offset
       local_triplets = total_values/this%process_grid%total_processors
       local_offset = local_triplets * (this%process_grid%global_rank)
       header_size = 4 * bytes_per_int
       IF (this%process_grid%global_rank .EQ. &
            & this%process_grid%total_processors - 1) THEN
          local_triplets = INT(total_values) - INT(local_offset)
       END IF
       local_offset = local_offset*(bytes_per_int*2+bytes_per_data) + &
            & header_size

       !! Do The Actual Reading
       CALL MPI_File_set_view(mpi_file_handler,local_offset,triplet_mpi_type,&
            & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
       IF (this%is_complex) THEN
          CALL ConstructTripletList(triplet_list_c, local_triplets)
          CALL MPI_File_read_all(mpi_file_handler, triplet_list_c%data, &
               & local_triplets, triplet_mpi_type, mpi_status,ierr)
       ELSE
          CALL ConstructTripletList(triplet_list_r, local_triplets)
          CALL MPI_File_read_all(mpi_file_handler, triplet_list_r%data, &
               & local_triplets, triplet_mpi_type, mpi_status,ierr)
       END IF
       CALL MPI_File_close(mpi_file_handler,ierr)
       CALL StopTimer("MPI Read Binary")

       IF (this%is_complex) THEN
          CALL FillMatrixFromTripletList(this,triplet_list_c)
          CALL DestructTripletList(triplet_list_c)
       ELSE
          CALL FillMatrixFromTripletList(this,triplet_list_r)
          CALL DestructTripletList(triplet_list_r)
       END IF
    END IF

  END SUBROUTINE ConstructMatrixFromBinary_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Save a distributed sparse matrix to a binary file.
  !> Faster than text, so this is good for check pointing.
  SUBROUTINE WriteMatrixToBinary_ps(this,file_name)
    !! Parameters
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! Local Data
    INTEGER :: triplet_mpi_type

    IF (this%is_complex) THEN
       triplet_mpi_type = GetMPITripletType_c()
       CALL WriteMatrixToBinary_psc(this, file_name, triplet_mpi_type)
    ELSE
       triplet_mpi_type = GetMPITripletType_r()
       CALL WriteMatrixToBinary_psr(this, file_name, triplet_mpi_type)
    END IF
  END SUBROUTINE WriteMatrixToBinary_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Implementation of write to binary.
  SUBROUTINE WriteMatrixToBinary_psr(this, file_name, triplet_mpi_type)
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !> The triplet type, which distinguishes real and complex triplets.
    INTEGER, INTENT(IN) :: triplet_mpi_type
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(Matrix_lsr) :: merged_local_data

    INCLUDE "distributed_includes/WriteMatrixToBinary.f90"

  END SUBROUTINE WriteMatrixToBinary_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Implementation of write to binary.
  SUBROUTINE WriteMatrixToBinary_psc(this, file_name, triplet_mpi_type)
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !> The triplet type, which distinguishes real and complex triplets.
    INTEGER, INTENT(IN) :: triplet_mpi_type
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(Matrix_lsc) :: merged_local_data

    INCLUDE "distributed_includes/WriteMatrixToBinary.f90"

  END SUBROUTINE WriteMatrixToBinary_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write a distributed sparse matrix to a matrix market file.
  !> Read \cite boisvert1996matrix for the details.
  SUBROUTINE WriteMatrixToMatrixMarket_ps(this,file_name)
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name

    IF (this%is_complex) THEN
       CALL WriteMatrixToMatrixMarket_psc(this, file_name)
    ELSE
       CALL WriteMatrixToMatrixMarket_psr(this, file_name)
    END IF
  END SUBROUTINE WriteMatrixToMatrixMarket_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write to matrix market implementation for real data.
  SUBROUTINE WriteMatrixToMatrixMarket_psr(this,file_name)
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(Matrix_lsr) :: merged_local_data

#include "distributed_includes/WriteToMatrixMarket.f90"

  END SUBROUTINE WriteMatrixToMatrixMarket_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Write to matrix market implementation for complex data.
  SUBROUTINE WriteMatrixToMatrixMarket_psc(this,file_name)
    !> The Matrix to write.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The name of the file to write to.
    CHARACTER(len=*), INTENT(IN) :: file_name
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(Matrix_lsc) :: merged_local_data

#define ISCOMPLEX
#include "distributed_includes/WriteToMatrixMarket.f90"
#undef ISCOMPLEX

  END SUBROUTINE WriteMatrixToMatrixMarket_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !> should pass in triplet lists with global coordinates. It doesn't matter
  !> where each triplet is stored, as long as global coordinates are given.
  SUBROUTINE FillMatrixFromTripletList_psr(this,triplet_list,preduplicated_in)
    !> The matrix to fill.
    TYPE(Matrix_ps) :: this
    !> The triplet list of values.
    TYPE(TripletList_r) :: triplet_list
    !> If lists are preduplicated across slices set this to true.
    LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    !! Local Data
    TYPE(Matrix_ps) :: temp_matrix
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix
    TYPE(Matrix_lsr) :: gathered_matrix
    !! Local Data
    TYPE(Permutation_t) :: basic_permutation
    TYPE(ReduceHelper_t) :: gather_helper
    REAL(NTREAL), PARAMETER :: threshold = 0.0
    LOGICAL :: preduplicated

    IF (this%is_complex) THEN
       CALL ConvertMatrixToReal(this, temp_matrix)
       CALL CopyMatrix(temp_matrix, this)
       CALL DestructMatrix(temp_matrix)
    END IF

    INCLUDE "distributed_includes/FillMatrixFromTripletList.f90"
  END SUBROUTINE FillMatrixFromTripletList_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine fills in a matrix based on local triplet lists. Each process
  !> should pass in triplet lists with global coordinates. It doesn't matter
  !> where each triplet is stored, as long as global coordinates are given.
  SUBROUTINE FillMatrixFromTripletList_psc(this,triplet_list,preduplicated_in)
    !> The matrix to fill.
    TYPE(Matrix_ps) :: this
    !> The triplet list of values.
    TYPE(TripletList_c) :: triplet_list
    !> If lists are preduplicated across slices set this to true.
    LOGICAL, INTENT(IN), OPTIONAL :: preduplicated_in
    !! Local Data
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix
    TYPE(Matrix_lsc) :: gathered_matrix
    !! Local Data
    TYPE(Matrix_ps) :: temp_matrix
    TYPE(Permutation_t) :: basic_permutation
    TYPE(ReduceHelper_t) :: gather_helper
    REAL(NTREAL), PARAMETER :: threshold = 0.0
    LOGICAL :: preduplicated

    IF (.NOT. this%is_complex) THEN
       CALL ConvertMatrixToComplex(this, temp_matrix)
       CALL CopyMatrix(temp_matrix, this)
       CALL DestructMatrix(temp_matrix)
    END IF

    INCLUDE "distributed_includes/FillMatrixFromTripletList.f90"
  END SUBROUTINE FillMatrixFromTripletList_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  SUBROUTINE FillMatrixIdentity_ps(this)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this

    IF (this%is_complex) THEN
       CALL FillMatrixIdentity_psc(this)
    ELSE
       CALL FillMatrixIdentity_psr(this)
    END IF

  END SUBROUTINE FillMatrixIdentity_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  SUBROUTINE FillMatrixIdentity_psr(this)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix

    INCLUDE "distributed_includes/FillMatrixIdentity.f90"

  END SUBROUTINE FillMatrixIdentity_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with the identity matrix.
  SUBROUTINE FillMatrixIdentity_psc(this)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: unsorted_triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix

    INCLUDE "distributed_includes/FillMatrixIdentity.f90"

  END SUBROUTINE FillMatrixIdentity_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill in the values of a distributed matrix with a permutation.
  !> If you don't specify permuterows, will default to permuting rows.
  SUBROUTINE FillMatrixPermutation_ps(this, permutation_vector, permute_rows_in)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Describes for each row/column, where it goes.
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    !> If true permute rows, false permute columns.
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
  !> Fill permutation implementation.
  SUBROUTINE FillMatrixPermutation_psr(this, permutation_vector, rows)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Describes for each row/column, where it goes.
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    !> If true permute rows, false permute columns.
    LOGICAL, INTENT(IN) :: rows
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: unsorted_triplet_list
    TYPE(TripletList_r) :: sorted_triplet_list
    TYPE(Matrix_lsr) :: local_matrix

    INCLUDE "distributed_includes/FillMatrixPermutation.f90"

  END SUBROUTINE FillMatrixPermutation_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill permutation implementation.
  SUBROUTINE FillMatrixPermutation_psc(this, permutation_vector, rows)
    !> The matrix being filled.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Describes for each row/column, where it goes.
    INTEGER, DIMENSION(:), INTENT(IN) :: permutation_vector
    !> If true permute rows, false permute columns.
    LOGICAL, INTENT(IN) :: rows
    !! Local Data
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: unsorted_triplet_list
    TYPE(TripletList_c) :: sorted_triplet_list
    TYPE(Matrix_lsc) :: local_matrix

    INCLUDE "distributed_includes/FillMatrixPermutation.f90"

  END SUBROUTINE FillMatrixPermutation_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !> Data is returned with absolute coordinates.
  SUBROUTINE GetMatrixTripletList_psr(this, triplet_list)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The list to fill.
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !! Local Data
    TYPE(Matrix_ps) :: working_matrix
    TYPE(Matrix_lsr) :: merged_local_data

    IF (this%is_complex) THEN
       CALL ConvertMatrixToReal(this, working_matrix)
    ELSE
       CALL CopyMatrix(this, working_matrix)
    END IF

    INCLUDE "distributed_includes/GetMatrixTripletList.f90"
  END SUBROUTINE GetMatrixTripletList_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extracts a triplet list of the data that is stored on this process.
  !> Data is returned with absolute coordinates.
  SUBROUTINE GetMatrixTripletList_psc(this, triplet_list)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The list to fill.
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !! Local Data
    TYPE(Matrix_ps) :: working_matrix
    TYPE(Matrix_lsc) :: merged_local_data

    IF (.NOT. this%is_complex) THEN
       CALL ConvertMatrixToComplex(this, working_matrix)
    ELSE
       CALL CopyMatrix(this, working_matrix)
    END IF

    INCLUDE "distributed_includes/GetMatrixTripletList.f90"
  END SUBROUTINE GetMatrixTripletList_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list. Block is
  !> defined by the row/column start/end values.
  !> This is slower than GetMatrixTripletList, because communication is required
  !> Data is returned with absolute coordinates.
  SUBROUTINE GetMatrixBlock_psr(this, triplet_list, start_row, end_row, &
       & start_column, end_column)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The list to fill.
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !> The starting row for data to store on this process.
    INTEGER :: start_row
    !> The ending row for data to store on this process.
    INTEGER :: end_row
    !> The starting col for data to store on this process
    INTEGER :: start_column
    !> The ending col for data to store on this process
    INTEGER :: end_column
    !! Local Data
    TYPE(Matrix_ps) :: working_matrix
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(TripletList_r) :: local_triplet_list
    !! Send Buffer
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    !! Receive Buffer
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! Temp Values
    TYPE(Triplet_r) :: temp_triplet
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_end_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_end_list
    !! Send Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
    !! Receive Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
    !! Temporary
    INTEGER :: II, PP
    INTEGER :: ierr

    IF (this%is_complex) THEN
       CALL ConvertMatrixToReal(this, working_matrix)
    ELSE
       CALL CopyMatrix(this, working_matrix)
    END IF

#define MPIDATATYPE MPINTREAL
#include "distributed_includes/GetMatrixBlock.f90"
#undef MPIDATATYPE

  END SUBROUTINE GetMatrixBlock_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract an arbitrary block of a matrix into a triplet list. Block is
  !> defined by the row/column start/end values.
  !> This is slower than GetMatrixTripletList, because communication is required
  !> Data is returned with absolute coordinates.
  SUBROUTINE GetMatrixBlock_psc(this, triplet_list, start_row, end_row, &
       & start_column, end_column)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The list to fill.
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !> The starting row for data to store on this process.
    INTEGER :: start_row
    !> The ending row for data to store on this process.
    INTEGER :: end_row
    !> The starting col for data to store on this process
    INTEGER :: start_column
    !> The ending col for data to store on this process
    INTEGER :: end_column
    !! Local Data
    TYPE(Matrix_ps) :: working_matrix
    TYPE(Matrix_lsc) :: merged_local_data
    TYPE(TripletList_c) :: local_triplet_list
    !! Send Buffer
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    !! Receive Buffer
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! Temp Values
    TYPE(Triplet_c) :: temp_triplet
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_start_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: row_end_list
    INTEGER, DIMENSION(:), ALLOCATABLE :: column_end_list
    !! Send Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
    !! Receive Buffer
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_proc
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
    !! Temporary
    INTEGER :: II, PP
    INTEGER :: ierr

    IF (.NOT. this%is_complex) THEN
       CALL ConvertMatrixToComplex(this, working_matrix)
    ELSE
       CALL CopyMatrix(this, working_matrix)
    END IF

#define MPIDATATYPE MPINTCOMPLEX
#include "distributed_includes/GetMatrixBlock.f90"
#undef MPIDATATYPE

  END SUBROUTINE GetMatrixBlock_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy an arbitrary slice from a matrix into a new smaller matrix.
  !> NTPoly only works with square matrices, so if the number of rows and
  !> columns is different the matrix is resized to the maximum size.
  SUBROUTINE GetMatrixSlice(this, submatrix, start_row, end_row, &
       & start_column, end_column)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The slice to fill.
    TYPE(Matrix_ps), INTENT(INOUT) :: submatrix
    !> The starting row to include in this matrix.
    INTEGER :: start_row
    !> The ending row to include in this matrix.
    INTEGER :: end_row
    !> The starting column to include in this matrix.
    INTEGER :: start_column
    !> The last column to include in this matrix.
    INTEGER :: end_column

    !! Get a triplet list with the values
    IF (this%is_complex) THEN
       CALL GetMatrixSlice_psc(this, submatrix, start_row, end_row, &
            & start_column, end_column)
    ELSE
       CALL GetMatrixSlice_psr(this, submatrix, start_row, end_row, &
            & start_column, end_column)
    END IF

  END SUBROUTINE GetMatrixSlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Implements slice matrix for real types.
  SUBROUTINE GetMatrixSlice_psr(this, submatrix, start_row, end_row, &
       & start_column, end_column)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The slice to fill.
    TYPE(Matrix_ps), INTENT(INOUT) :: submatrix
    !> The starting row to include in this matrix.
    INTEGER :: start_row
    !> The ending row to include in this matrix.
    INTEGER :: end_row
    !> The starting column to include in this matrix.
    INTEGER :: start_column
    !> The last column to include in this matrix.
    INTEGER :: end_column

#define TLISTTYPE TripletList_r
#define TTYPE Triplet_r
#include "distributed_includes/SliceMatrix.f90"
#undef TLISTTYPE
#undef TTYPE

  END SUBROUTINE GetMatrixSlice_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Implements slice matrix for complex types.
  SUBROUTINE GetMatrixSlice_psc(this, submatrix, start_row, end_row, &
       & start_column, end_column)
    !> The distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The slice to fill.
    TYPE(Matrix_ps), INTENT(INOUT) :: submatrix
    !> The starting row to include in this matrix.
    INTEGER :: start_row
    !> The ending row to include in this matrix.
    INTEGER :: end_row
    !> The starting column to include in this matrix.
    INTEGER :: start_column
    !> The last column to include in this matrix.
    INTEGER :: end_column

#define TLISTTYPE TripletList_c
#define TTYPE Triplet_c
#include "distributed_includes/SliceMatrix.f90"
#undef TLISTTYPE
#undef TTYPE

  END SUBROUTINE GetMatrixSlice_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the actual dimension of the matrix.
  PURE FUNCTION GetMatrixActualDimension_ps(this) RESULT(DIMENSION)
    !! Parameters
    !> The matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Dimension of the matrix
    INTEGER :: DIMENSION
    DIMENSION = this%actual_matrix_dimension
  END FUNCTION GetMatrixActualDimension_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the logical dimension of the matrix.
  !> Includes padding.
  PURE FUNCTION GetMatrixLogicalDimension_ps(this) RESULT(DIMENSION)
    !> The matrix.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Dimension of the matrix
    INTEGER :: DIMENSION
    DIMENSION = this%logical_matrix_dimension
  END FUNCTION GetMatrixLogicalDimension_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out information about a distributed sparse matrix.
  !> Sparsity, and load balancing information.
  SUBROUTINE PrintMatrixInformation_ps(this)
    !> This the matrix to print information about.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !! Local Data
    INTEGER :: min_size, max_size
    REAL(NTREAL) :: sparsity

    CALL GetMatrixLoadBalance(this,min_size,max_size)
    sparsity = REAL(GetMatrixSize(this),KIND=NTREAL) / &
         & (REAL(this%actual_matrix_dimension,KIND=NTREAL)**2)

    CALL WriteHeader("Load_Balance")
    CALL EnterSubLog
    CALL WriteListElement(key="min_size", value=min_size)
    CALL WriteListElement(key="max_size", value=max_size)
    CALL ExitSubLog
    CALL WriteElement(key="Dimension",value=this%actual_matrix_dimension)
    CALL WriteElement(key="Sparsity", value=sparsity)
  END SUBROUTINE PrintMatrixInformation_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a distributed sparse matrix.
  !> This is a serial print routine, and should probably only be used for debug
  !> purposes.
  SUBROUTINE PrintMatrix_ps(this, file_name_in)
    !> The matrix to print.
    TYPE(Matrix_ps) :: this
    !> Optionally, you can pass a file to print to instead of the console.
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
  !> Print matrix implementation (real).
  SUBROUTINE PrintMatrix_psr(this, file_name_in)
    !> The matrix to print.
    TYPE(Matrix_ps) :: this
    !> Optionally, you can pass a file to print to instead of the console.
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Temporary Variables
    TYPE(Matrix_lsr) :: local_mat

    INCLUDE "distributed_includes/PrintMatrix.f90"
  END SUBROUTINE PrintMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print matrix implementation (complex).
  SUBROUTINE PrintMatrix_psc(this, file_name_in)
    !> The matrix to print.
    TYPE(Matrix_ps) :: this
    !> Optionally, you can pass a file to print to instead of the console.
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Temporary Variables
    TYPE(Matrix_lsc) :: local_mat

    INCLUDE "distributed_includes/PrintMatrix.f90"
  END SUBROUTINE PrintMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A utility routine that filters a sparse matrix.
  !> All (absolute) values below the threshold are set to zero.
  SUBROUTINE FilterMatrix_ps(this, threshold)
    !> The matrix to filter.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Threshold (absolute) values below this are filtered
    REAL(NTREAL), INTENT(IN) :: threshold

    IF (this%is_complex) THEN
       CALL FilterMatrix_psc(this, threshold)
    ELSE
       CALL FilterMatrix_psr(this, threshold)
    END IF
  END SUBROUTINE FilterMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Filter matrix implementation (real).
  SUBROUTINE FilterMatrix_psr(this, threshold)
    !> The matrix to filter.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Threshold (absolute) values below this are filtered
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary

    INCLUDE "distributed_includes/FilterMatrix.f90"
  END SUBROUTINE FilterMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Filter matrix implementation (real).
  SUBROUTINE FilterMatrix_psc(this, threshold)
    !> The matrix to filter.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Threshold (absolute) values below this are filtered
    REAL(NTREAL), INTENT(IN) :: threshold
    !! Local Variables
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: new_list
    TYPE(Triplet_c) :: temporary

    INCLUDE "distributed_includes/FilterMatrix.f90"
  END SUBROUTINE FilterMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the total number of non-zero entries in the distributed sparse matrix.
  FUNCTION GetMatrixSize_ps(this) RESULT(total_size)
    !> The matrix to calculate the number of non-zero entries of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The number of non-zero entries in the matrix.
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

    total_size = INT(temp_size,kind=c_long)

  END FUNCTION GetMatrixSize_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get a measure of how load balanced this matrix is. For each process, the
  !> number of non-zero entries is calculated. Then, this function returns
  !> the max and min of those values.
  SUBROUTINE GetMatrixLoadBalance_ps(this, min_size, max_size)
    !> The matrix to compute the measure on.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The minimum entries contained on a single process.
    INTEGER, INTENT(OUT) :: min_size
    !> The maximum entries contained on a single process.
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
    CALL MPI_Allreduce(local_size,max_size,1,MPINTINTEGER,MPI_MAX,&
         & this%process_grid%within_slice_comm, ierr)
    CALL MPI_Allreduce(local_size,min_size,1,MPINTINTEGER,MPI_MIN,&
         & this%process_grid%within_slice_comm, ierr)

  END SUBROUTINE GetMatrixLoadBalance_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix. Note that this is a pure transpose, there is
  !> no complex conjugate performed.
  SUBROUTINE TransposeMatrix_ps(AMat, TransMat)
    !> The matrix to transpose.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> TransMat = A^T .
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat

    IF (AMat%is_complex) THEN
       CALL TransposeMatrix_psc(AMat, TransMat)
    ELSE
       CALL TransposeMatrix_psr(AMat, TransMat)
    END IF

  END SUBROUTINE TransposeMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose implementation (real).
  SUBROUTINE TransposeMatrix_psr(AMat, TransMat)
    !> The matrix to transpose.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> TransMat = A^T .
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat
    !! Local Variables
    TYPE(TripletList_r) :: triplet_list
    TYPE(TripletList_r) :: new_list
    TYPE(Triplet_r) :: temporary, temporary_t

    INCLUDE "distributed_includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose implementation (complex).
  SUBROUTINE TransposeMatrix_psc(AMat, TransMat)
    !> The matrix to transpose.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> TransMat = A^T .
    TYPE(Matrix_ps), INTENT(OUT) :: TransMat
    !! Local Variables
    TYPE(TripletList_c) :: triplet_list
    TYPE(TripletList_c) :: new_list
    TYPE(Triplet_c) :: temporary, temporary_t

    INCLUDE "distributed_includes/TransposeMatrix.f90"

  END SUBROUTINE TransposeMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Every value in the matrix is changed into its complex conjugate.
  PURE SUBROUTINE ConjugateMatrix_ps(this)
    !> The matrix to compute the complex conjugate of.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !! Local Variables
    TYPE(Matrix_lsc) :: local_matrix

    IF (this%is_complex) THEN
       CALL MergeMatrixLocalBlocks(this, local_matrix)
       CALL ConjugateMatrix(local_matrix)
       CALL SplitMatrixToLocalBlocks(this, local_matrix)
       CALL DestructMatrix(local_matrix)
    END IF

  END SUBROUTINE ConjugateMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split the current communicator, and give each group a complete copy of this
  SUBROUTINE CommSplitMatrix_ps(this, split_mat, my_color, split_slice)
    !> The matrix to split.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A copy of the matrix hosted on a small process grid.
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    !> Distinguishes between the two groups.
    INTEGER, INTENT(OUT) :: my_color
    !> If we split along the slice direction, this is True
    LOGICAL, INTENT(OUT) :: split_slice

    IF (this%is_complex) THEN
       CALL CommSplitMatrix_psc(this, split_mat, my_color, split_slice)
    ELSE
       CALL CommSplitMatrix_psr(this, split_mat, my_color, split_slice)
    END IF

  END SUBROUTINE CommSplitMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split implementation for real data.
  SUBROUTINE CommSplitMatrix_psr(this, split_mat, my_color, split_slice)
    !> The matrix to split.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A copy of the matrix hosted on a small process grid.
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    !> Distinguishes between the two groups.
    INTEGER, INTENT(OUT) :: my_color
    !> If we split along the slice direction, this is True.
    LOGICAL, INTENT(OUT) :: split_slice
    !! For Data Redistribution
    TYPE(TripletList_r) :: full_list, new_list
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: send_list

    INCLUDE "distributed_includes/CommSplitMatrix.f90"

  END SUBROUTINE CommSplitMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split implementation for complex data.
  SUBROUTINE CommSplitMatrix_psc(this, split_mat, my_color, split_slice)
    !> The matrix to split.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A copy of the matrix hosted on a small process grid.
    TYPE(Matrix_ps), INTENT(INOUT) :: split_mat
    !> Distinguishes between the two groups.
    INTEGER, INTENT(OUT) :: my_color
    !> If we split along the slice direction, this is True.
    LOGICAL, INTENT(OUT) :: split_slice
    !! For Data Redistribution
    TYPE(TripletList_c) :: full_list, new_list
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: send_list

    INCLUDE "distributed_includes/CommSplitMatrix.f90"

  END SUBROUTINE CommSplitMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute the data in a matrix based on row, column list
  !> This will redistribute the data so that the local data are entries in
  !> the rows and columns list. The order of the row list and column list matter
  !> because local data is filled in the same order.
  SUBROUTINE RedistributeData_psr(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !> The matrix to redistribute
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Lookup describing how data is distributed.
    INTEGER, DIMENSION(:), INTENT(IN) :: index_lookup
    !> Reverse Lookup describing how data is distributed.
    INTEGER, DIMENSION(:), INTENT(IN) :: reverse_index_lookup
    !> The current triplet list of global coordinates.
    TYPE(TripletList_r), INTENT(IN) :: initial_triplet_list
    !> returns an allocated triplet list with local coordinates in sorted order.
    TYPE(TripletList_r), INTENT(OUT) :: sorted_triplet_list
    !! Local Data
    TYPE(TripletList_r) :: gathered_list
    TYPE(TripletList_r), DIMENSION(this%process_grid%slice_size) :: &
         & send_triplet_lists
    TYPE(Triplet_r) :: temp_triplet

    INCLUDE "distributed_includes/RedistributeData.f90"

  END SUBROUTINE RedistributeData_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute the data in a matrix based on row, column list
  !> This will redistribute the data so that the local data are entries in
  !> the rows and columns list. The order of the row list and column list matter
  !> because local data is filled in the same order.
  SUBROUTINE RedistributeData_psc(this,index_lookup,reverse_index_lookup,&
       & initial_triplet_list,sorted_triplet_list)
    !> The matrix to redistribute
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> Lookup describing how data is distributed.
    INTEGER, DIMENSION(:), INTENT(IN) :: index_lookup
    !> Reverse Lookup describing how data is distributed.
    INTEGER, DIMENSION(:), INTENT(IN) :: reverse_index_lookup
    !> The current triplet list of global coordinates.
    TYPE(TripletList_c), INTENT(IN) :: initial_triplet_list
    !> returns an allocated triplet list with local coordinates in sorted order.
    TYPE(TripletList_c), INTENT(OUT) :: sorted_triplet_list
    !! Local Data
    TYPE(TripletList_c) :: gathered_list
    TYPE(TripletList_c), DIMENSION(this%process_grid%slice_size) :: &
         & send_triplet_lists
    TYPE(Triplet_c) :: temp_triplet

    INCLUDE "distributed_includes/RedistributeData.f90"

  END SUBROUTINE RedistributeData_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Calculate a matrix size that can be divided by the number of processors.
  PURE FUNCTION CalculateScaledDimension(this, matrix_dim) RESULT(scaled_dim)
    !> The matrix we are calculating for.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The dimension of the actual matrix.
    INTEGER, INTENT(IN) :: matrix_dim
    !> A new dimension which includes padding.
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
  PURE SUBROUTINE SplitMatrixToLocalBlocks_psr(this, matrix_to_split)
    !> The distributed sparse matrix to split into.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The matrix to split up.
    TYPE(Matrix_lsr), INTENT(IN) :: matrix_to_split

#define LOCALMATRIX this%local_data_r
#include "distributed_includes/SplitMatrixToLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE SplitMatrixToLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a local matrix, and use it to fill the local block matrix structure.
  PURE SUBROUTINE SplitMatrixToLocalBlocks_psc(this, matrix_to_split)
    !> The distributed sparse matrix to split into.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The matrix to split up.
    TYPE(Matrix_lsc), INTENT(IN) :: matrix_to_split

#define LOCALMATRIX this%local_data_c
#include "distributed_includes/SplitMatrixToLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE SplitMatrixToLocalBlocks_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge together the local matrix blocks into one big matrix.
  PURE SUBROUTINE MergeMatrixLocalBlocks_psr(this, merged_matrix)
    !> The distributed sparse matrix to merge from.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The merged matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: merged_matrix

#define LOCALMATRIX this%local_data_r
#include "distributed_includes/MergeMatrixLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE MergeMatrixLocalBlocks_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Merge together the local matrix blocks into one big matrix.
  !! @param[inout] this the distributed sparse matrix.
  !! @param[inout] merged_matrix the merged matrix.
  PURE SUBROUTINE MergeMatrixLocalBlocks_psc(this, merged_matrix)
    !> The distributed sparse matrix to merge from.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The merged matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: merged_matrix

#define LOCALMATRIX this%local_data_c
#include "distributed_includes/MergeMatrixLocalBlocks.f90"
#undef LOCALMATRIX

  END SUBROUTINE MergeMatrixLocalBlocks_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the current matrix to a real type matrix.
  SUBROUTINE ConvertMatrixToReal(in, out)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: in
    !> Real version of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: out
    LOGICAL, PARAMETER :: convert_to_complex = .FALSE.
    !! Local Variables
    TYPE(Matrix_lsc) :: local_matrix
    TYPE(Matrix_lsr) :: converted_matrix

    INCLUDE "distributed_includes/ConvertMatrixType.f90"
  END SUBROUTINE ConvertMatrixToReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Converts the current matrix to a complex type matrix.
  SUBROUTINE ConvertMatrixToComplex(in, out)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: in
    !> Complex version of the matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: out
    LOGICAL, PARAMETER :: convert_to_complex = .TRUE.
    !! Local Variables
    TYPE(Matrix_lsr) :: local_matrix
    TYPE(Matrix_lsc) :: converted_matrix

    INCLUDE "distributed_includes/ConvertMatrixType.f90"
  END SUBROUTINE ConvertMatrixToComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Change the size of a matrix.
  !! If the new size is smaller, then values outside that range are deleted.
  !! IF the new size is bigger, zero padding is applied.
  !! Warning: this requires a full data redistribution.
  SUBROUTINE ResizeMatrix(this, new_size)
    !> The matrix to resize.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The new size of the matrix.
    INTEGER, INTENT(IN) :: new_size

    IF (this%is_complex) THEN
       CALL ResizeMatrix_psc(this, new_size)
    ELSE
       CALL ResizeMatrix_psr(this, new_size)
    END IF
  END SUBROUTINE ResizeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Change the size of a matrix implementation (real).
  SUBROUTINE ResizeMatrix_psr(this, new_size)
    !> The matrix to resize.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The new size of the matrix.
    INTEGER, INTENT(IN) :: new_size
    !! Local Variables
    TYPE(TripletList_r) :: tlist, pruned
    TYPE(Triplet_r) :: temp

    INCLUDE "distributed_includes/ResizeMatrix.f90"
  END SUBROUTINE ResizeMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Change the size of a matrix implementation (real).
  SUBROUTINE ResizeMatrix_psc(this, new_size)
    !> The matrix to resize.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> The new size of the matrix.
    INTEGER, INTENT(IN) :: new_size
    !! Local Variables
    TYPE(TripletList_c) :: tlist, pruned
    TYPE(Triplet_c) :: temp

    INCLUDE "distributed_includes/ResizeMatrix.f90"
  END SUBROUTINE ResizeMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This subroutine gathers the entire matrix into a local matrix on the
  !! given process. This routine is used when printing, but also is useful for
  !! debugging.
  SUBROUTINE GatherMatrixToProcess_psr(this, local_mat, proc_id)
    !> The matrix to gather.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The full matrix, stored in a local matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: local_mat
    !> Which process to gather on.
    INTEGER, INTENT(IN) :: proc_id
    !! Local Variables
    TYPE(TripletList_r) :: tlist, sorted
    TYPE(TripletList_r), DIMENSION(:), ALLOCATABLE :: slist

    INCLUDE "distributed_includes/GatherMatrixToProcess.f90"
  END SUBROUTINE GatherMatrixToProcess_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This subroutine gathers the entire matrix into a local matrix on the
  !! given process. This routine is used when printing, but also is useful for
  !! debugging.
  SUBROUTINE GatherMatrixToProcess_psc(this, local_mat, proc_id)
    !> The matrix to gather.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The full matrix, stored in a local matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: local_mat
    !> Which process to gather on.
    INTEGER, INTENT(IN) :: proc_id
    !! Local Variables
    TYPE(TripletList_c) :: tlist, sorted
    TYPE(TripletList_c), DIMENSION(:), ALLOCATABLE :: slist

    INCLUDE "distributed_includes/GatherMatrixToProcess.f90"
  END SUBROUTINE GatherMatrixToProcess_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixModule
