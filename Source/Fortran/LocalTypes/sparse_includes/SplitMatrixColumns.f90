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
     CALL ConstructEmptyMatrix(split_list(split_counter), this%rows, lcolumns)
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
