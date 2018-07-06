  !! Local Data
  INTEGER, DIMENSION(:), ALLOCATABLE :: row_lookup
  INTEGER, DIMENSION(:), ALLOCATABLE :: column_lookup
  INTEGER, DIMENSION(:), ALLOCATABLE :: location_list_within_slice
  !! Temporary Values
  INTEGER :: row_size, column_size
  INTEGER :: temp_row, temp_column
  INTEGER :: process_id
  INTEGER :: counter

  CALL StartTimer("Redistribute")

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
     CALL ConstructTripletList(send_triplet_lists(counter))
  END DO
  DO counter = 1, initial_triplet_list%CurrentSize
     process_id = location_list_within_slice(counter)
     CALL GetTripletAt(initial_triplet_list, counter, temp_triplet)
     CALL AppendToTripletList(send_triplet_lists(process_id+1), temp_triplet)
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
  CALL StartTimer("SortTripletList")
  CALL SortTripletList(gathered_list, this%local_columns, this%local_rows, &
       & sorted_triplet_list)
  CALL StopTimer("SortTripletList")

  !! Cleanup
  DO counter = 1, this%process_grid%slice_size
     CALL DestructTripletList(send_triplet_lists(counter))
  END DO
  DEALLOCATE(row_lookup)
  DEALLOCATE(column_lookup)
  DEALLOCATE(location_list_within_slice)
  CALL DestructTripletList(gathered_list)

  CALL StopTimer("Redistribute")
