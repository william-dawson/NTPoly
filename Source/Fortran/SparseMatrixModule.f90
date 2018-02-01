!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling locally stored CSR matrices.
MODULE SparseMatrixModule
  USE DataTypesModule, ONLY : NTREAL
  USE MatrixMarketModule
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule, ONLY: TripletList_t, SortTripletList, &
       & ConstructTripletList, DestructTripletList, SetTripletAt, &
       & AppendToTripletList, SymmetrizeTripletList
  USE TripletModule, ONLY : Triplet_t
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a CSR matrix.
  TYPE, PUBLIC :: SparseMatrix_t
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index !< Outer indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: inner_index !< Inner indices
     REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: values !< Values
     INTEGER :: rows !< Matrix dimension: rows
     INTEGER :: columns !< Matrix dimension: columns
  END TYPE SparseMatrix_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Construct/Destruct
  PUBLIC :: ConstructEmptySparseMatrix
  PUBLIC :: ConstructSparseMatrixFromFile
  PUBLIC :: ConstructZeroSparseMatrix
  PUBLIC :: ConstructFromTripletList
  PUBLIC :: DestructSparseMatrix
  PUBLIC :: CopySparseMatrix
  !! Basic Accessors
  PUBLIC :: GetRows
  PUBLIC :: GetColumns
  PUBLIC :: ExtractRow
  PUBLIC :: ExtractColumn
  !! Helper routines
  PUBLIC :: SplitSparseMatrixColumns
  PUBLIC :: ComposeSparseMatrixColumns
  PUBLIC :: TransposeSparseMatrix
  PUBLIC :: PrintSparseMatrix
  PUBLIC :: MatrixToTripletList
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmptySparseMatrix(this, columns, rows)
    !! Parameters
    TYPE(SparseMatrix_t),INTENT(OUT)  :: this
    INTEGER, INTENT(IN) :: columns, rows

    CALL DestructSparseMatrix(this)
    this%rows = rows
    this%columns = columns
    ALLOCATE(this%outer_index(this%columns+1))
    this%outer_index = 0
  END SUBROUTINE ConstructEmptySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  SUBROUTINE ConstructSparseMatrixFromFile(this, file_name)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(OUT) :: this
    CHARACTER(len=*), INTENT(IN)   :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    INTEGER :: temp_rows, temp_columns, temp_total_values
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    CHARACTER(len=81) :: input_buffer
    INTEGER :: file_handler
    TYPE(Triplet_t) :: temporary
    INTEGER :: counter
    LOGICAL :: found_comment_line
    LOGICAL :: error_occured
    file_handler = 16

    OPEN(file_handler,file=file_name,status='old')

    !! Parse the header.
    READ(file_handler,fmt='(A)') input_buffer
    error_occured = ParseMMHeader(input_buffer, sparsity_type, data_type, &
         & pattern_type)

    !! Extra Comment Lines
    found_comment_line = .TRUE.
    DO WHILE(found_comment_line)
       !read(file_handler,*) input_buffer
       READ(file_handler,fmt='(A)') input_buffer
       IF (.NOT. input_buffer(1:1) .EQ. '%') THEN
          found_comment_line = .FALSE.
       END IF
    END DO

    !! Main data
    READ(input_buffer,*) temp_rows, temp_columns, temp_total_values
    CALL ConstructTripletList(triplet_list)

    !! Read Values
    DO counter = 1, temp_total_values
       READ(file_handler,*) temporary%index_row, temporary%index_column, &
            & temporary%point_value
       CALL AppendToTripletList(triplet_list,temporary)
    END DO
    CLOSE(file_handler)

    CALL SymmetrizeTripletList(triplet_list, pattern_type)
    CALL SortTripletList(triplet_list,temp_columns,sorted_triplet_list)
    CALL ConstructFromTripletList(this, sorted_triplet_list, temp_rows, &
         & temp_columns)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE ConstructSparseMatrixFromFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix with zero values in it.
  !! @param[out] this the matrix being constructed
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructZeroSparseMatrix(this,rows,columns)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(OUT) :: this
    INTEGER, INTENT(IN)             :: rows, columns

    !! Allocate
    CALL ConstructEmptySparseMatrix(this,columns,rows)
    ALLOCATE(this%inner_index(0))
    ALLOCATE(this%values(0))

  END SUBROUTINE ConstructZeroSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructFromTripletList(this,triplet_list,rows,columns)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(OUT) :: this
    TYPE(TripletList_t), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN)             :: rows, columns
    INTEGER :: outer_array_ptr
    INTEGER :: values_counter

    IF (ALLOCATED(this%outer_index)) THEN
       CALL DestructSparseMatrix(this)
    END IF

    this%rows = rows
    this%columns = columns

    !! Allocate
    ALLOCATE(this%outer_index(this%columns+1))
    this%outer_index = 0
    ALLOCATE(this%inner_index(triplet_list%CurrentSize))
    ALLOCATE(this%values(triplet_list%CurrentSize))

    !! Insert Values
    outer_array_ptr = 1
    DO values_counter = 1, triplet_list%CurrentSize
       !! Moving on to the next column?
       DO WHILE(.NOT. triplet_list%data(values_counter)%index_column .EQ. &
            & outer_array_ptr)
          outer_array_ptr = outer_array_ptr + 1
          this%outer_index(outer_array_ptr+1) = this%outer_index(outer_array_ptr)
       END DO
       this%outer_index(outer_array_ptr+1)=this%outer_index(outer_array_ptr+1)+1
       !! Insert inner index and value
       this%inner_index(values_counter) = &
            & triplet_list%data(values_counter)%index_row
       this%values(values_counter) = &
            & triplet_list%data(values_counter)%point_value
    END DO

    !! Fill In The Rest Of The Outer Values
    DO outer_array_ptr = outer_array_ptr+2, this%columns+1
       this%outer_index(outer_array_ptr) = this%outer_index(outer_array_ptr-1)
    END DO

  END SUBROUTINE ConstructFromTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE DestructSparseMatrix(this)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%outer_index)) THEN
       DEALLOCATE(this%outer_index)
    END IF
    IF (ALLOCATED(this%inner_index)) THEN
       DEALLOCATE(this%inner_index)
    END IF
    IF (ALLOCATED(this%values)) THEN
       DEALLOCATE(this%values)
    END IF
  END SUBROUTINE DestructSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopySparseMatrix(matA, matB)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: matA
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matB

    CALL DestructSparseMatrix(matB)
    matB = matA
  END SUBROUTINE CopySparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetRows(this) RESULT(rows)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    INTEGER :: rows
    rows = this%rows
  END FUNCTION GetRows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetColumns(this) RESULT(columns)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    INTEGER :: columns
    columns = this%columns
  END FUNCTION GetColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractRow(this, row_number, row_out)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    TYPE(SparseMatrix_t), INTENT(INOUT) :: row_out
    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: counter
    TYPE(SparseMatrix_t) :: temp, temp_c

    CALL TransposeSparseMatrix(this,temp)
    CALL ExtractColumn(temp, row_number, temp_c)
    CALL TransposeSparseMatrix(temp_c,row_out)
    CALL DestructSparseMatrix(temp_c)
    CALL DestructSparseMatrix(temp)
  END SUBROUTINE ExtractRow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractColumn(this, column_number, column_out)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    TYPE(SparseMatrix_t), INTENT(INOUT) :: column_out
    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: counter

    !! Allocate Memory
    CALL ConstructEmptySparseMatrix(column_out, 1, this%rows)
    start_index = this%outer_index(column_number)
    number_of_values = this%outer_index(column_number+1) - &
         & this%outer_index(column_number)
    ALLOCATE(column_out%inner_index(number_of_values))
    ALLOCATE(column_out%values(number_of_values))

    !! Copy Values
    column_out%outer_index(1) = 0
    column_out%outer_index(2) = number_of_values
    DO counter=1, number_of_values
      column_out%inner_index(counter) = this%inner_index(start_index+counter)
      column_out%values(counter) = this%values(start_index+counter)
    END DO
  END SUBROUTINE ExtractColumn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE TransposeSparseMatrix(this, matT)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN)  :: this
    TYPE(SparseMatrix_t), INTENT(INOUT) :: matT
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    TYPE(TripletList_t) :: sorted_triplet_list
    !! Temporary variables and counters
    INTEGER :: counter
    INTEGER :: temp_int
    INTEGER :: temp_columns, temp_rows
    INTEGER :: size_of_this

    CALL DestructSparseMatrix(matT)

    temp_rows = this%columns
    temp_columns = this%rows

    !! Construct a transposed triplet list
    size_of_this = this%outer_index(this%columns+1)
    CALL ConstructTripletList(triplet_list,size_of_this)

    CALL MatrixToTripletList(this,triplet_list)
    DO counter = 1,triplet_list%CurrentSize
       temp_int = triplet_list%data(counter)%index_column
       triplet_list%data(counter)%index_column = &
            & triplet_list%data(counter)%index_row
       triplet_list%data(counter)%index_row = temp_int
    END DO

    !! Build a new matrix from the triplet list
    CALL SortTripletList(triplet_list,temp_columns,sorted_triplet_list, &
         & bubble_in=.FALSE.)
    CALL ConstructFromTripletList(matT, sorted_triplet_list,&
         & temp_rows,temp_columns)

    !! Cleanup
    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE TransposeSparseMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !! the first matrix are followed by the columns of the matrices in the list.
  !! @param[in] mat_list list of matrices to compose.
  !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
  PURE SUBROUTINE ComposeSparseMatrixColumns(mat_list, out_matrix)
    !! Parameters
    TYPE(SparseMatrix_t), DIMENSION(:), INTENT(IN) :: mat_list
    TYPE(SparseMatrix_t), INTENT(INOUT) :: out_matrix
    !! Local Variables
    INTEGER :: total_columns, total_values
    INTEGER :: inner_start, inner_length
    INTEGER :: outer_start, outer_length
    INTEGER :: outer_offset
    INTEGER :: counter
    INTEGER :: size_of_mat

    CALL DestructSparseMatrix(out_matrix)

    !! Figure Out The Sizes
    total_columns = 0
    total_values  = 0
    DO counter = LBOUND(mat_list,dim=1), UBOUND(mat_list,dim=1)
       total_columns = total_columns + mat_list(counter)%columns
       size_of_mat = mat_list(counter)%outer_index(mat_list(counter)%columns+1)
       total_values  = total_values + size_of_mat
    END DO

    !! Allocate The Space
    CALL ConstructEmptySparseMatrix(out_matrix,total_columns,mat_list(1)%rows)
    ALLOCATE(out_matrix%inner_index(total_values))
    ALLOCATE(out_matrix%values(total_values))

    !! Fill In The Values
    inner_start = 1
    outer_start = 1
    outer_offset = 0
    DO counter = LBOUND(mat_list,dim=1),UBOUND(mat_list,dim=1)
       !! Inner indices and values
       size_of_mat = mat_list(counter)%outer_index(mat_list(counter)%columns+1)
       inner_length = size_of_mat
       out_matrix%inner_index(inner_start:inner_start+inner_length-1) = &
            & mat_list(counter)%inner_index
       out_matrix%values(inner_start:inner_start+inner_length-1) = &
            & mat_list(counter)%values
       inner_start = inner_start + inner_length
       !! Outer Indices
       outer_length = mat_list(counter)%columns+1
       out_matrix%outer_index(outer_start:outer_start+outer_length-1) = &
            & mat_list(counter)%outer_index + outer_offset
       outer_start = outer_start + outer_length - 1
       outer_offset = out_matrix%outer_index(outer_start)
    END DO

  END SUBROUTINE ComposeSparseMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Take a matrix, and split into into small blocks.
  !! @param[in] this matrix to perform this operation on.
  !! @param[in] num_blocks number of blocks to split into.
  !! @param[out] split_list 1D array of blocks.
  !! @param[out] block_offsets_out the offsets used for splitting.
  PURE SUBROUTINE SplitSparseMatrixColumns(this, num_blocks, split_list, &
       & block_offsets_out)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: num_blocks
    TYPE(SparseMatrix_t), DIMENSION(num_blocks), INTENT(OUT) :: split_list
    INTEGER, DIMENSION(num_blocks+1), INTENT(OUT), OPTIONAL :: block_offsets_out
    !! Local Data
    INTEGER, DIMENSION(num_blocks) :: block_sizes
    INTEGER, DIMENSION(num_blocks+1) :: block_offsets
    INTEGER :: split_divisor
    !! Counters
    INTEGER :: split_counter
    !! Temporary variables
    INTEGER :: loffset, lcolumns, linner_offset, total_values

    !! Handle trivial case
    IF (num_blocks .EQ. 1) THEN
       CALL CopySparseMatrix(this,split_list(1))
       block_offsets(1) = 1
       block_offsets(2) = this%columns+1
    ELSE
       !! Compute the sizes of each block
       split_divisor = this%columns/num_blocks
       block_sizes = split_divisor
       !! Handle an uneven split
       block_sizes(num_blocks) = this%columns - split_divisor*(num_blocks-1)
       !! Offsets
       block_offsets(1) = 1
       DO split_counter = 2, num_blocks+1
          block_offsets(split_counter) = block_offsets(split_counter-1) + &
               & block_sizes(split_counter-1)
       END DO

       !! Split up the columns
       DO split_counter = 1, num_blocks
          !! Temporary variables
          loffset = block_offsets(split_counter)
          lcolumns = block_sizes(split_counter)
          linner_offset = this%outer_index(loffset)+1
          !! Construct
          CALL ConstructEmptySparseMatrix(split_list(split_counter), &
               & columns=lcolumns, rows=this%rows)
          !! Copy Outer Index
          split_list(split_counter)%outer_index =        &
               & this%outer_index(loffset:loffset+lcolumns)
          split_list(split_counter)%outer_index =        &
               & split_list(split_counter)%outer_index -    &
               & split_list(split_counter)%outer_index(1)
          total_values = split_list(split_counter)%outer_index(lcolumns+1)
          !! Copy Inner Indices and Values
          IF (total_values .GT. 0) THEN
             ALLOCATE(split_list(split_counter)%inner_index(total_values))
             split_list(split_counter)%inner_index = &
                  & this%inner_index(linner_offset:linner_offset+total_values-1)
             ALLOCATE(split_list(split_counter)%values(total_values))
             split_list(split_counter)%values = &
                  & this%values(linner_offset:linner_offset+total_values-1)
          END IF
       END DO
    END IF
    IF (PRESENT(block_offsets_out)) THEN
       block_offsets_out = block_offsets
    END IF
  END SUBROUTINE SplitSparseMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  !! @param[in] this the matrix to construct the triplet list from.
  !! @param[out] triplet_list the triplet list we created.
  PURE SUBROUTINE MatrixToTripletList(this, triplet_list)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    TYPE(TripletList_t), INTENT(INOUT) :: triplet_list
    !! Helper variables
    INTEGER :: outer_counter, inner_counter
    INTEGER :: elements_per_inner
    INTEGER :: total_counter
    TYPE(Triplet_t) :: temporary
    INTEGER :: size_of_this

    size_of_this = this%outer_index(this%columns+1)
    CALL ConstructTripletList(triplet_list,size_of_this)

    total_counter = 1
    DO outer_counter = 1, this%columns
       elements_per_inner = this%outer_index(outer_counter+1) - &
            & this%outer_index(outer_counter)
       DO inner_counter = 1, elements_per_inner
          temporary%index_column = outer_counter
          temporary%index_row = this%inner_index(total_counter)
          temporary%point_value = this%values(total_counter)
          triplet_list%data(total_counter) = temporary
          total_counter = total_counter + 1
       END DO
    END DO
  END SUBROUTINE MatrixToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a sparse matrix.
  !! @param[in] this the matrix to be printed.
  !! @param[in] file_name_in optionally, you can pass a file to print to.
  SUBROUTINE PrintSparseMatrix(this, file_name_in)
    !! Parameters
    TYPE(SparseMatrix_t), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TripletList_t) :: triplet_list
    INTEGER :: file_handler
    INTEGER :: counter
    INTEGER :: size_of_this

    !! Process Optional Parameters
    IF (PRESENT(file_name_in)) THEN
       file_handler = 16
       OPEN(unit = file_handler, file = file_name_in)
    ELSE
       file_handler = 6
    END IF

    !! Print
    CALL MatrixToTripletList(this,triplet_list)

    size_of_this = this%outer_index(this%columns+1)

    WRITE(file_handler,'(A)') "%%MatrixMarket matrix coordinate real general"
    WRITE(file_handler,'(A)') "%"
    WRITE(file_handler,*) this%rows, this%columns, size_of_this
    DO counter = 1,size_of_this
       WRITE(file_handler,*) triplet_list%data(counter)%index_row, &
            & triplet_list%data(counter)%index_column, &
            & triplet_list%data(counter)%point_value
    END DO

    IF (PRESENT(file_name_in)) THEN
       CLOSE(file_handler)
    END IF
    CALL DestructTripletList(triplet_list)
  END SUBROUTINE PrintSparseMatrix
END MODULE SparseMatrixModule
