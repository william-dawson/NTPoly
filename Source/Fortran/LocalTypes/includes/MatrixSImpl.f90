!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Internal only. Create a sparse matrix with a certain number of columns
  !! and rows. Will allocate storage for the outer values, nothing else.
  !! @param[out] this the matrix being created. It will have the outer
  !! index allocated, but nothing else.
  !! @param[in] columns number of matrix columns.
  !! @param[in] rows number of matrix rows.
  PURE SUBROUTINE ConstructEmptyMatrix(this, columns, rows)
    !! Parameters
    TYPE(SMTYPE),INTENT(OUT)  :: this
    INTEGER, INTENT(IN) :: columns, rows

    CALL DestructMatrix(this)
    this%rows = rows
    this%columns = columns
    ALLOCATE(this%outer_index(this%columns+1))
    this%outer_index = 0
  END SUBROUTINE ConstructEmptyMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a sparse matrix by reading in a matrix market file.
  !! @param[out] this the matrix being constructed.
  !! @param[in] file_name name of the file.
  SUBROUTINE ConstructMatrixFromFile(this, file_name)
    !! Parameters
    TYPE(SMTYPE), INTENT(OUT) :: this
    CHARACTER(len=*), INTENT(IN)   :: file_name
    !! About the matrix market file.
    INTEGER :: sparsity_type, data_type, pattern_type
    !! Local Data
    INTEGER :: temp_rows, temp_columns, temp_total_values
    TYPE(TLISTTYPE) :: triplet_list
    TYPE(TLISTTYPE) :: sorted_triplet_list
    CHARACTER(len=81) :: input_buffer
    INTEGER :: file_handler
    TYPE(TTYPE) :: temporary
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
    CALL SortTripletList(triplet_list, temp_columns, temp_rows, &
         & sorted_triplet_list)
    CALL ConstructMatrixFromTripletList(this, sorted_triplet_list, temp_rows, &
         & temp_columns)

    CALL DestructTripletList(triplet_list)
    CALL DestructTripletList(sorted_triplet_list)
  END SUBROUTINE ConstructMatrixFromFile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix with zero values in it.
  !! @param[out] this the matrix being constructed
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructZeroMatrix(this,rows,columns)
    !! Parameters
    TYPE(SMTYPE), INTENT(OUT) :: this
    INTEGER, INTENT(IN)             :: rows, columns

    !! Allocate
    CALL ConstructEmptyMatrix(this,columns,rows)
    ALLOCATE(this%inner_index(0))
    ALLOCATE(this%values(0))

  END SUBROUTINE ConstructZeroMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a sparse matrix from a \b SORTED triplet list.
  !! The triplet list must be sorted to efficiently fill in the matrix. This
  !! constructor assumes \b you have already sorted the triplet list.
  !! @param[out] this the matrix being constructed
  !! @param[in] triplet_list a list of triplet values. They must be sorted.
  !! @param[in] rows number of matrix rows
  !! @param[in] columns number of matrix columns
  PURE SUBROUTINE ConstructMatrixFromTripletList(this,triplet_list,rows,columns)
    !! Parameters
    TYPE(SMTYPE), INTENT(OUT) :: this
    TYPE(TLISTTYPE), INTENT(IN) :: triplet_list
    INTEGER, INTENT(IN)             :: rows, columns
    INTEGER :: outer_array_ptr
    INTEGER :: values_counter

    IF (ALLOCATED(this%outer_index)) THEN
       CALL DestructMatrix(this)
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

  END SUBROUTINE ConstructMatrixFromTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Explicitly destruct a sparse matrix.
  !! @param[inout] this the matrix to free up
  PURE SUBROUTINE DestructMatrix(this)
    !! Parameters
    TYPE(SMTYPE), INTENT(INOUT) :: this

    IF (ALLOCATED(this%outer_index)) THEN
       DEALLOCATE(this%outer_index)
    END IF
    IF (ALLOCATED(this%inner_index)) THEN
       DEALLOCATE(this%inner_index)
    END IF
    IF (ALLOCATED(this%values)) THEN
       DEALLOCATE(this%values)
    END IF
  END SUBROUTINE DestructMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy a sparse matrix in a safe way.
  !! @param[in] matA matrix to copy
  !! @param[inout] matB = matA
  PURE SUBROUTINE CopyMatrix(matA, matB)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: matA
    TYPE(SMTYPE), INTENT(INOUT) :: matB

    CALL DestructMatrix(matB)
    matB = matA
  END SUBROUTINE CopyMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of rows of a matrix.
  !! @param[in] this the matrix.
  !! @result number of rows.
  PURE FUNCTION GetMatrixRows(this) RESULT(rows)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER :: rows
    rows = this%rows
  END FUNCTION GetMatrixRows
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of columns of a matrix.
  !! @param[in] this the matrix.
  !! @result number of columns.
  PURE FUNCTION GetMatrixColumns(this) RESULT(columns)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER :: columns
    columns = this%columns
  END FUNCTION GetMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a row from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] row_number the row to extract
  !! @param[out] row_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixRow(this, row_number, row_out)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: row_number
    TYPE(SMTYPE), INTENT(INOUT) :: row_out
    !! Temporary Variables
    INTEGER :: values_found
    DATATYPE, DIMENSION(:), ALLOCATABLE :: value_buffer
    INTEGER :: total_counter, elements_per_inner
    INTEGER :: outer_counter
    INTEGER :: inner_counter

    !! Fill a value buffer
    CALL ConstructEmptyMatrix(row_out,this%columns,1)
    ALLOCATE(value_buffer(this%columns))
    values_found = 0
    total_counter = 1
    row_out%outer_index(1) = 0
    DO outer_counter = 1, this%columns
       row_out%outer_index(outer_counter+1) = &
            & row_out%outer_index(outer_counter+1) + &
            & row_out%outer_index(outer_counter)
       elements_per_inner = this%outer_index(outer_counter+1) - &
            & this%outer_index(outer_counter)
       DO inner_counter = 1, elements_per_inner
          IF (this%inner_index(total_counter) .EQ. row_number) THEN
             values_found = values_found + 1
             value_buffer(values_found) = this%values(total_counter)
             row_out%outer_index(outer_counter+1) = &
                  & row_out%outer_index(outer_counter+1) + 1
          END IF
          total_counter = total_counter + 1
       END DO
    END DO

    !! Copy To Actual Matrix
    ALLOCATE(row_out%inner_index(values_found))
    row_out%inner_index = 1
    ALLOCATE(row_out%values(values_found))
    row_out%values = value_buffer(:values_found)

    !! Cleanup
    DEALLOCATE(value_buffer)
  END SUBROUTINE ExtractMatrixRow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Extract a column from the matrix into the compressed vector representation.
  !! @param[in] this the matrix to extract from.
  !! @param[in] column_number the row to extract
  !! @param[out] column_out the matrix representing that row
  PURE SUBROUTINE ExtractMatrixColumn(this, column_number, column_out)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: column_number
    TYPE(SMTYPE), INTENT(INOUT) :: column_out
    !! Local variables
    INTEGER :: number_of_values
    INTEGER :: start_index
    INTEGER :: counter

    !! Allocate Memory
    CALL ConstructEmptyMatrix(column_out, 1, this%rows)
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
  END SUBROUTINE ExtractMatrixColumn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a sparse matrix and return it in a separate matrix.
  !! The current implementation has you go from matrix to triplet list,
  !! triplet list to transposed triplet list. The triplet list must then be
  !! sorted and then the return matrix is constructed.
  !! @param[in] this the matrix to be transposed.
  !! @param[out] matT the input matrix transposed.
  PURE SUBROUTINE TransposeMatrix(this, matT)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN)  :: this
    TYPE(SMTYPE), INTENT(INOUT) :: matT
    !! Local Data
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    !! Temporary Variables
    INTEGER :: II, JJ
    INTEGER :: inner_index, insert_pt, this_offset
    INTEGER :: num_values, elements_per_inner

    !! Allocate New Matrix
    num_values = this%outer_index(this%columns+1)
    CALL ConstructEmptyMatrix(matT,this%rows,this%columns)
    ALLOCATE(matT%inner_index(num_values))
    ALLOCATE(matT%values(num_values))

    !! Temporary Arrays
    ALLOCATE(values_per_row(this%rows))
    ALLOCATE(offset_array(this%rows))

    !! Count the values per row
    values_per_row = 0
    DO II = 1, num_values
       inner_index = this%inner_index(II)
       values_per_row(inner_index) = values_per_row(inner_index) + 1
    END DO
    offset_array(1) = 0
    DO II = 2, this%rows
       offset_array(II) = offset_array(II-1) + values_per_row(II-1)
    END DO

    !! Insert
    matT%outer_index(:this%rows) = offset_array(:this%rows)
    matT%outer_index(this%rows+1) = offset_array(this%rows) + &
         & values_per_row(this%rows)
    DO II = 1, this%columns
       elements_per_inner = this%outer_index(II+1) - this%outer_index(II)
       this_offset = this%outer_index(II)
       DO JJ = 1, elements_per_inner
          inner_index = this%inner_index(this_offset+JJ)
          insert_pt = offset_array(inner_index)+1
          matT%inner_index(insert_pt) = II
          matT%values(insert_pt) = this%values(this_offset+JJ)
          offset_array(inner_index) = offset_array(inner_index) +1
       END DO
    END DO

    !! Cleanup
    DEALLOCATE(values_per_row)
    DEALLOCATE(offset_array)
  END SUBROUTINE TransposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !! to another.
  !! @param[in] mat_array 2d array of matrices to compose.
  !! @param[in] block_rows the number of rows of the array of blocks.
  !! @param[in] block_columns the number of columns of the array of blocks.
  !! @param[out] out_matrix the composed matrix.
  PURE SUBROUTINE ComposeMatrix(mat_array, block_rows, block_columns, &
       & out_matrix)
    !! Parameters
    TYPE(SMTYPE), DIMENSION(block_rows,block_columns), INTENT(IN) :: &
         & mat_array
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(SMTYPE), INTENT(INOUT) :: out_matrix
    !! Local Data
    TYPE(SMTYPE), DIMENSION(block_columns) :: merged_columns
    TYPE(SMTYPE) :: Temp
    TYPE(SMTYPE), DIMENSION(block_rows,block_columns) :: mat_t
    INTEGER :: II, JJ

    !! First transpose the matrices
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL TransposeMatrix(mat_array(II,JJ), mat_t(II,JJ))
       END DO
    END DO

    !! Next merge the columns
    DO JJ = 1, block_columns
       CALL ComposeMatrixColumns(mat_t(:,JJ), Temp)
       CALL TransposeMatrix(Temp, merged_columns(JJ))
    END DO

    !! Final Merge
    CALL ComposeMatrixColumns(merged_columns, out_matrix)

    !! Cleanup
    DO JJ = 1, block_columns
       DO II = 1, block_rows
          CALL DestructMatrix(mat_t(II,JJ))
       END DO
    END DO
    DO JJ = 1, block_columns
       CALL DestructMatrix(merged_columns(JJ))
    END DO

  END SUBROUTINE ComposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big Matrix C = [Matrix 1 | Matrix 1, ...] where the columns of
  !! the first matrix are followed by the columns of the matrices in the list.
  !! @param[in] mat_list list of matrices to compose.
  !! @param[out] out_matrix = [Matrix 1 | Matrix 2, ...] .
  PURE SUBROUTINE ComposeMatrixColumns(mat_list, out_matrix)
    !! Parameters
    TYPE(SMTYPE), DIMENSION(:), INTENT(IN) :: mat_list
    TYPE(SMTYPE), INTENT(INOUT) :: out_matrix
    !! Local Variables
    INTEGER :: total_columns, total_values
    INTEGER :: inner_start, inner_length
    INTEGER :: outer_start, outer_length
    INTEGER :: outer_offset
    INTEGER :: counter
    INTEGER :: size_of_mat

    CALL DestructMatrix(out_matrix)

    !! Figure Out The Sizes
    total_columns = 0
    total_values  = 0
    DO counter = LBOUND(mat_list,dim=1), UBOUND(mat_list,dim=1)
       total_columns = total_columns + mat_list(counter)%columns
       size_of_mat = mat_list(counter)%outer_index(mat_list(counter)%columns+1)
       total_values  = total_values + size_of_mat
    END DO

    !! Allocate The Space
    CALL ConstructEmptyMatrix(out_matrix,total_columns,mat_list(1)%rows)
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

  END SUBROUTINE ComposeMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  !! @param[in] this the matrix to split.
  !! @param[in] block_rows number of rows to split the matrix into.
  !! @param[in] block_columns number of columns to split the matrix into.
  !! @param[out] split_array a COLUMNxROW array for the output to go into.
  !! @param[in] block_size_row_in specifies the size of the  rows.
  !! @param[in] block_size_column_in specifies the size of the columns.
  SUBROUTINE SplitMatrix(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: block_rows, block_columns
    TYPE(SMTYPE), DIMENSION(:,:), INTENT(INOUT) :: split_array
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in
    !! Local Data
    INTEGER, DIMENSION(block_rows) :: block_size_row
    INTEGER, DIMENSION(block_columns) :: block_size_column
    !! Temporary Variables
    TYPE(SMTYPE), DIMENSION(block_columns) :: column_split
    TYPE(SMTYPE), DIMENSION(block_rows) :: row_split
    TYPE(SMTYPE) :: Temp
    INTEGER :: divisor_row, divisor_column
    INTEGER :: II, JJ

    !! Calculate the split sizes
    IF (PRESENT(block_size_row_in)) THEN
       block_size_row = block_size_row_in
    ELSE
       divisor_row = this%rows/block_rows
       block_size_row = divisor_row
       block_size_row(block_rows) = this%rows - divisor_row*(block_rows-1)
    END IF
    IF (PRESENT(block_size_column_in)) THEN
       block_size_column = block_size_column_in
    ELSE
       divisor_column = this%columns/block_columns
       block_size_column = divisor_column
       block_size_column(block_columns) = this%columns - &
            & divisor_column*(block_columns-1)
    END IF

    !! First split by columns which is easy with the CSR format
    CALL StartTimer("Split Column")
    CALL SplitMatrixColumns(this, block_columns, block_size_column, &
         & column_split)
    CALL StopTimer("Split Column")

    !! Now Split By Rows
    CALL StartTimer("Split Row")
    DO JJ = 1, block_columns
       CALL TransposeMatrix(column_split(JJ), Temp)
       CALL SplitMatrixColumns(Temp, block_rows, block_size_row, &
            & row_split)
       !! Copy into output array
       DO II = 1, block_rows
          CALL TransposeMatrix(row_split(II), split_array(II,JJ))
       END DO
    END DO
    CALL StopTimer("Split Row")

    !! Cleanup
    CALL DestructMatrix(Temp)
    DO II = 1, block_rows
       CALL DestructMatrix(row_split(II))
    END DO
    DO II = 1, block_columns
       CALL DestructMatrix(column_split(II))
    END DO

  END SUBROUTINE SplitMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a matrix into into small blocks based on the specified offsets.
  !! @param[in] this matrix to perform this operation on.
  !! @param[in] num_blocks number of blocks to split into
  !! @param[out] block_sizes the sizes used for splitting.
  !! @param[out] split_list 1D array of blocks.
  PURE SUBROUTINE SplitMatrixColumns(this, num_blocks, block_sizes, &
       & split_list)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: num_blocks
    INTEGER, DIMENSION(num_blocks), INTENT(IN) :: block_sizes
    TYPE(SMTYPE), DIMENSION(num_blocks), INTENT(INOUT) :: split_list
    !! Local Data
    INTEGER, DIMENSION(num_blocks+1) :: block_offsets
    !! Counters
    INTEGER :: split_counter
    !! Temporary variables
    INTEGER :: loffset, lcolumns, linner_offset, total_values

    !! Compute Offsets
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
       CALL ConstructEmptyMatrix(split_list(split_counter), &
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
  END SUBROUTINE SplitMatrixColumns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list from a matrix.
  !! @param[in] this the matrix to construct the triplet list from.
  !! @param[out] triplet_list the triplet list we created.
  PURE SUBROUTINE MatrixToTripletList(this, triplet_list)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    TYPE(TLISTTYPE), INTENT(INOUT) :: triplet_list
    !! Helper variables
    INTEGER :: outer_counter, inner_counter
    INTEGER :: elements_per_inner
    INTEGER :: total_counter
    TYPE(TTYPE) :: temporary
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
  SUBROUTINE PrintMatrix(this, file_name_in)
    !! Parameters
    TYPE(SMTYPE), INTENT(IN) :: this
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: file_name_in
    !! Local Data
    TYPE(TLISTTYPE) :: triplet_list
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
  END SUBROUTINE PrintMatrix
