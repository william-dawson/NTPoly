  !> Construct a triplet list.
  !! @param[inout] this the triplet list to construct.
  !! @param[in] size_in the length of the triplet list (optional, default=0).
  PURE FUNCTION ConstructTripletList(size_in) RESULT(this)
    !! Parameters
    TYPE(TLISTTYPE) :: this
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in
    !! Local data
    INTEGER :: size

    IF (PRESENT(size_in)) THEN
       size = size_in
    ELSE
       size = 0
    END IF

    IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
    this%CurrentSize = size

    ALLOCATE(this%data(size))
  END FUNCTION ConstructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  !! @param[inout] this the triplet list to destruct.
  PURE SUBROUTINE DestructTripletList(this)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: this

    IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
    this%CurrentSize = 0
  END SUBROUTINE DestructTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  !! @param[inout] this the triplet list to resize.
  !! @param[in] size to resize to.
  PURE SUBROUTINE ResizeTripletList(this, size)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(TTYPE), DIMENSION(:), ALLOCATABLE :: temporary_data

    !! Temporary copy
    ALLOCATE(temporary_data(this%CurrentSize))
    temporary_data = this%DATA(:this%CurrentSize)

    !! Create new memory
    IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
    ALLOCATE(this%DATA(size))

    !! Copy back
    this%DATA(:this%CurrentSize) = temporary_data

    !! Cleanup
    DEALLOCATE(temporary_data)
  END SUBROUTINE ResizeTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  !! @param[inout] this the triplet list to append to.
  !! @param[in] triplet_value the value to append.
  PURE SUBROUTINE AppendToTripletList(this, triplet_value)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: this
    TYPE(TTYPE), INTENT(IN)        :: triplet_value
    !! Local data
    INTEGER :: new_size

    !! First, check if we need to allocate more memory
    IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
       IF (SIZE(this%data) .EQ. 0) THEN
          new_size = 1
       ELSE IF (SIZE(this%data) .EQ. 1) THEN
          new_size = 2
       ELSE
          new_size = INT(SIZE(this%data)*1.5)
       END IF
       CALL ResizeTripletList(this,new_size)
    END IF

    !! Append
    this%CurrentSize = this%CurrentSize+1
    this%DATA(this%CurrentSize) = triplet_value

  END SUBROUTINE AppendToTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> (Just for a related project)
  PURE SUBROUTINE AccumulateTripletList(this, triplet_value)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: this
    TYPE(TTYPE), INTENT(IN)        :: triplet_value
    !! Local data
    INTEGER :: new_size
    INTEGER :: idata

    DO idata = 1, this%CurrentSize
       IF ((this%DATA(idata)%index_row == triplet_value%index_row) .AND. &
            &    (this%DATA(idata)%index_column == triplet_value%index_column)) THEN
          this%DATA(idata)%point_value = this%DATA(idata)%point_value + triplet_value%point_value
          RETURN
       END IF
    END DO

    !! First, check if we need to allocate more memory
    IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
       IF (SIZE(this%data) .EQ. 0) THEN
          new_size = 1
       ELSE IF (SIZE(this%data) .EQ. 1) THEN
          new_size = 2
       ELSE
          new_size = INT(SIZE(this%data)*1.5)
       END IF
       CALL ResizeTripletList(this,new_size)
    END IF

    !! Append
    this%CurrentSize = this%CurrentSize+1
    this%DATA(this%CurrentSize) = triplet_value

  END SUBROUTINE AccumulateTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  !! @param[inout] this the triplet list to set.
  !! @param[in] index the index at which to set the triplet.
  !! @param[in] triplet_value the value of the triplet to set.
  PURE SUBROUTINE SetTripletAt(this,index,triplet_value)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: this
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    TYPE(TTYPE), INTENT(IN)        :: triplet_value

    this%data(index) = triplet_value
  END SUBROUTINE SetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  !! @param[in] this the triplet list to get the value from.
  !! @param[in] index the index from which to get the triplet.
  !! @param[out] triplet_value the extracted triplet value.
  PURE SUBROUTINE GetTripletAt(this,index,triplet_value)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(IN) :: this
    INTEGER(kind=c_int), INTENT(IN) :: index
    TYPE(TTYPE), INTENT(OUT)    :: triplet_value

    triplet_value = this%data(index)
  END SUBROUTINE GetTripletAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !! Implementation is based on bucket sort. This is why it needs the number of
  !! matrix columns. Bubble sort is used within a bucket.
  !! @param[in] input_list list to be sorted.
  !! @param[in] matrix_columns this is the highest column value in the list.
  !! @param[in] matrix_rows this is the highest row value in the list.
  !! @param[in] bubble_in false if you don't need the final bubble sort.
  !! @param[out] sorted_list a now sorted version of the list. This routine
  !! will allocate it.
  PURE SUBROUTINE SortTripletList(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TLISTTYPE), INTENT(OUT) :: sorted_list
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    LOGICAL :: bubble
    TYPE(TTYPE) :: temporary
    LOGICAL :: swap_occured
    INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
    INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
    !! Counters and temporary variables
    INTEGER :: counter
    INTEGER :: temp_index
    INTEGER :: alloc_stat
    INTEGER :: list_length

    IF (PRESENT(bubble_in)) THEN
       bubble = bubble_in
    ELSE
       bubble = .TRUE.
    END IF

    list_length = input_list%CurrentSize

    IF (bubble .AND. list_length .GT. matrix_rows*matrix_columns*0.1) THEN
       CALL SortDenseTripletList(input_list, matrix_columns, matrix_rows, &
            & sorted_list)
    ELSE
       !! Data Allocation
       sorted_list = ConstructTripletList(list_length)
       ALLOCATE(values_per_row(matrix_columns), stat=alloc_stat)
       ALLOCATE(offset_array(matrix_columns), stat=alloc_stat)
       ALLOCATE(inserted_per_row(matrix_columns), stat=alloc_stat)

       !! Initial one dimensional sort
       values_per_row = 0
       inserted_per_row = 0

       !! Do a first pass bucket sort
       DO counter = 1, input_list%CurrentSize
          values_per_row(input_list%data(counter)%index_column) = &
               & values_per_row(input_list%data(counter)%index_column) + 1
       END DO
       offset_array(1) = 1
       DO counter = 2, UBOUND(offset_array,dim=1)
          offset_array(counter) = offset_array(counter-1) + &
               & values_per_row(counter-1)
       END DO
       DO counter = 1, input_list%CurrentSize
          temp_index = input_list%data(counter)%index_column
          sorted_list%data(offset_array(temp_index)+inserted_per_row(temp_index))=&
               & input_list%data(counter)
          inserted_per_row(temp_index) = inserted_per_row(temp_index) + 1
       END DO

       !! Finish with bubble sort
       !! Not necessary for transposing or unpacking.
       swap_occured = .TRUE.
       IF (bubble) THEN
          DO WHILE (swap_occured .EQV. .TRUE.)
             swap_occured = .FALSE.
             DO counter = 2, sorted_list%CurrentSize
                IF (CompareTriplets(sorted_list%data(counter-1), &
                     & sorted_list%data(counter))) THEN
                   temporary = sorted_list%data(counter)
                   sorted_list%data(counter) = sorted_list%data(counter-1)
                   sorted_list%data(counter-1) = temporary
                   swap_occured = .TRUE.
                END IF
             END DO
          END DO
       END IF

       !! Cleanup
       DEALLOCATE(values_per_row)
       DEALLOCATE(offset_array)
       DEALLOCATE(inserted_per_row)
    END IF

  END SUBROUTINE SortTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  !! @param[in] triplet_list list to get the size of.
  !! @return list_size the number of entries in the triplet list.
  PURE FUNCTION GetTripletListSize(triplet_list) RESULT(list_size)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(IN)  :: triplet_list
    INTEGER :: list_size
    list_size = triplet_list%CurrentSize
  END FUNCTION GetTripletListSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !! Takes in a list of triplet lists, one list for each processor. Then the
  !! all to all redistribution is performed along the given communicator.
  !! @param[in] triplet_lists a list of triplet lists, one for each process.
  !! @param[inout] comm the mpi communicator to redistribute along.
  !! @param[out] local_data_out the resulting local triplet list.
  SUBROUTINE RedistributeTripletLists(triplet_lists, comm, local_data_out)
    !! Parameters
    TYPE(TLISTTYPE), DIMENSION(:), INTENT(IN) :: triplet_lists
    INTEGER, INTENT(INOUT) :: comm
    TYPE(TLISTTYPE), INTENT(INOUT) :: local_data_out
    !! Local Data - Offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_process
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_offsets
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_process
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_offsets
    !! Local Data - Send/Recv Buffers
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
    DATATYPE, DIMENSION(:), ALLOCATABLE :: send_buffer_val
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
    INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
    DATATYPE, DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    !! ETC
    TYPE(TTYPE) :: temp_triplet
    INTEGER :: num_processes
    INTEGER :: counter, inner_counter, insert_pt
    INTEGER :: mpi_error

    !! Allocate Size Buffers
    CALL MPI_COMM_SIZE(comm, num_processes, mpi_error)
    ALLOCATE(send_per_process(num_processes))
    ALLOCATE(send_offsets(num_processes))
    ALLOCATE(recv_per_process(num_processes))
    ALLOCATE(recv_offsets(num_processes))

    !! Figure Out How Much Data Gets Sent
    DO counter = 1, num_processes
       send_per_process(counter) = triplet_lists(counter)%CurrentSize
    END DO
    send_offsets(1) = 0
    DO counter = 2, num_processes
       send_offsets(counter) = send_offsets(counter-1) + &
            & send_per_process(counter-1)
    END DO

    !! Figure Out How Much Data Gets Received
    CALL MPI_ALLTOALL(send_per_process, 1, MPI_INT, recv_per_process, 1, &
         & MPI_INT, comm, mpi_error)
    recv_offsets(1) = 0
    DO counter = 2, num_processes
       recv_offsets(counter) = recv_offsets(counter-1) + &
            & recv_per_process(counter-1)
    END DO

    !! Allocate And Fill Send Buffers
    ALLOCATE(send_buffer_row(SUM(send_per_process)))
    ALLOCATE(send_buffer_col(SUM(send_per_process)))
    ALLOCATE(send_buffer_val(SUM(send_per_process)))
    ALLOCATE(recv_buffer_row(SUM(recv_per_process)))
    ALLOCATE(recv_buffer_col(SUM(recv_per_process)))
    ALLOCATE(recv_buffer_val(SUM(recv_per_process)))

    !! Fill Send Buffer
    insert_pt = 1
    DO counter = 1, num_processes
       DO inner_counter = 1, triplet_lists(counter)%CurrentSize
          CALL GetTripletAt(triplet_lists(counter), inner_counter, temp_triplet)
          send_buffer_row(insert_pt) = temp_triplet%index_row
          send_buffer_col(insert_pt) = temp_triplet%index_column
          send_buffer_val(insert_pt) = temp_triplet%point_value
          insert_pt = insert_pt + 1
       END DO
    END DO

    !! Do Actual Send
    CALL StartTimer("AllToAllV")
    CALL MPI_Alltoallv(send_buffer_col, send_per_process, send_offsets, &
         & MPI_INT, recv_buffer_col, recv_per_process, recv_offsets, MPI_INT, &
         & comm, mpi_error)
    CALL MPI_Alltoallv(send_buffer_row, send_per_process, send_offsets, &
         & MPI_INT, recv_buffer_row, recv_per_process, recv_offsets, MPI_INT, &
         & comm, mpi_error)
    CALL MPI_Alltoallv(send_buffer_val, send_per_process, send_offsets, &
         & MPIDATATYPE, recv_buffer_val, recv_per_process, recv_offsets, &
         & MPIDATATYPE, comm, mpi_error)
    CALL StopTimer("AllToAllV")

    !! Unpack Into The Output Triplet List
    local_data_out = ConstructTripletList(size_in=SUM(recv_per_process))
    DO counter = 1, SUM(recv_per_process)
       local_data_out%data(counter)%index_column = recv_buffer_col(counter)
       local_data_out%data(counter)%index_row = recv_buffer_row(counter)
       local_data_out%data(counter)%point_value = recv_buffer_val(counter)
    END DO

    !! Cleanup
    DEALLOCATE(send_per_process)
    DEALLOCATE(send_offsets)
    DEALLOCATE(recv_per_process)
    DEALLOCATE(recv_offsets)
    DEALLOCATE(send_buffer_row)
    DEALLOCATE(send_buffer_col)
    DEALLOCATE(send_buffer_val)
    DEALLOCATE(recv_buffer_row)
    DEALLOCATE(recv_buffer_col)
    DEALLOCATE(recv_buffer_val)

  END SUBROUTINE RedistributeTripletLists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !! Frequently, we have a triplet list that comes from the global matrix which
  !! we would like to shift into a local matrix. In that case, just pass
  !! the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList(triplet_list, row_shift, column_shift)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(INOUT) :: triplet_list
    INTEGER, INTENT(IN) :: row_shift
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: counter

    !! Loop
    DO counter = 1, triplet_list%CurrentSize
       triplet_list%data(counter)%index_row = &
            triplet_list%data(counter)%index_row + row_shift
       triplet_list%data(counter)%index_column = &
            triplet_list%data(counter)%index_column + column_shift
    END DO

  END SUBROUTINE ShiftTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !! dense.
  !! @param[in] input_list the list
  !! @param[in] matrix_columns for the corresponding matrix.
  !! @param[in] matrix_rows for the corresponding matrix.
  !! @param[out] sorted_list sorted and ready to use for building matrices.
  PURE SUBROUTINE SortDenseTripletList(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !! Parameters
    TYPE(TLISTTYPE), INTENT(IN)  :: input_list
    INTEGER, INTENT(IN) :: matrix_columns
    INTEGER, INTENT(IN) :: matrix_rows
    TYPE(TLISTTYPE), INTENT(OUT) :: sorted_list
    !! Local Variables
    DATATYPE, DIMENSION(:,:), ALLOCATABLE :: value_buffer
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: dirty_buffer
    INTEGER :: list_length
    INTEGER :: row, col, ind
    INTEGER :: II, JJ

    !! Setup Memory
    ALLOCATE(value_buffer(matrix_rows,matrix_columns))
    ALLOCATE(dirty_buffer(matrix_rows,matrix_columns))
    value_buffer = 0
    dirty_buffer = 0
    list_length = input_list%CurrentSize
    sorted_list = ConstructTripletList(list_length)

    !! Unpack
    DO II = 1, list_length
       row = input_list%data(II)%index_row
       col = input_list%data(II)%index_column
       value_buffer(row,col) = input_list%data(II)%point_value
       dirty_buffer(row,col) = 1
    END DO

    !! Repack
    ind = 1
    DO JJ = 1, matrix_columns
       DO II = 1, matrix_rows
          IF (dirty_buffer(II,JJ) .EQ. 1) THEN
             sorted_list%data(ind)%index_row = II
             sorted_list%data(ind)%index_column = JJ
             sorted_list%data(ind)%point_value = value_buffer(II,JJ)
             ind = ind + 1
          END IF
       END DO
    END DO

    !! Cleanup
    DEALLOCATE(value_buffer)
    DEALLOCATE(dirty_buffer)

  END SUBROUTINE SortDenseTripletList
