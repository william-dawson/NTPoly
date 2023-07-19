  !! Local Data
  INTEGER, DIMENSION(num_blocks+1) :: block_offsets
  !! Temporary variables
  INTEGER :: II, JJ, loffset, lcolumns, linner_offset, total_values

  !! Compute Offsets
  block_offsets(1) = 1
  DO II = 2, num_blocks + 1
     block_offsets(II) = block_offsets(II - 1) + block_sizes(II - 1)
  END DO

  !! Split up the columns
  DO II = 1, num_blocks
     !! Temporary variables
     loffset = block_offsets(II)
     lcolumns = block_sizes(II)
     linner_offset = this%outer_index(loffset) + 1

     !! Construct
     CALL ConstructEmptyMatrix(split_list(II), this%rows, lcolumns)

     !! Copy Outer Index
     split_list(II)%outer_index(:) = &
          & this%outer_index(loffset:loffset + lcolumns)
     split_list(II)%outer_index(:) = split_list(II)%outer_index(:) - &
          & split_list(II)%outer_index(1)
     total_values = split_list(II)%outer_index(lcolumns + 1)

     !! Copy Inner Indices and Values
     IF (total_values .GT. 0) THEN
        ALLOCATE(split_list(II)%inner_index(total_values))
        ALLOCATE(split_list(II)%values(total_values))
        DO JJ = 1, total_values
          split_list(II)%inner_index(JJ) = this%inner_index(linner_offset + JJ - 1)
          split_list(II)%values(JJ) = this%values(linner_offset + JJ - 1)
        END DO
     END IF
  END DO
