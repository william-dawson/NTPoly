  !! Local MPI Variables
  INTEGER :: mpi_file_handler
  INTEGER :: message_status(MPI_STATUS_SIZE)
  INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
  !! Local Data
  INTEGER :: triplet_list_string_length
  INTEGER(KIND = MPI_OFFSET_KIND) :: header_size
  INTEGER(KIND = MPI_OFFSET_KIND) :: write_offset
  INTEGER(KIND = MPI_OFFSET_KIND) :: header_offset
  INTEGER(KIND = MPI_OFFSET_KIND), PARAMETER :: zero_size = 0
  !! Strings
  CHARACTER(LEN = :), ALLOCATABLE :: header_line1
  CHARACTER(LEN = :), ALLOCATABLE :: header_line2
  CHARACTER(LEN = :), ALLOCATABLE :: write_buffer
  !! Temporary Values
  INTEGER :: II, OFF_JJ
  INTEGER :: NEW_LINE_LENGTH
  CHARACTER(LEN = MAX_LINE_LENGTH*2) :: temp_string1
  CHARACTER(LEN = MAX_LINE_LENGTH) :: temp_string2
  CHARACTER(LEN = MAX_LINE_LENGTH) :: temp_string3
  INTEGER :: temp_length
  INTEGER :: bytes_per_character
  INTEGER :: ierr

  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data)

  CALL MPI_Type_size(MPI_CHARACTER, bytes_per_character, ierr)

  !! Create the matrix size line
  NEW_LINE_LENGTH = LEN(new_LINE('A'))
#ifdef ISCOMPLEX
  WRITE(temp_string1, '(A)') &
       & "%%MatrixMarket matrix coordinate complex general" &
       & // new_LINE('A') // "%" // new_LINE('A')
#else
  WRITE(temp_string1, '(A)') "%%MatrixMarket matrix coordinate real general" &
       & // new_LINE('A') // "%" // new_LINE('A')
#endif
  ALLOCATE(CHARACTER(LEN = LEN_TRIM(temp_string1)) :: header_line1)
  header_line1(:) = TRIM(temp_string1)

  CALL WriteMMSize(temp_string2, this%actual_matrix_dimension, &
       & this%actual_matrix_dimension, GetMatrixSize(this))
  ALLOCATE(CHARACTER(&
       & LEN = LEN_TRIM(temp_string2) + NEW_LINE_LENGTH + 1) :: header_line2)
  WRITE(header_line2,*) TRIM(temp_string2) // new_LINE('A')

  header_size = LEN(header_line1) + LEN(header_line2)

  !! Local Data
  CALL MatrixToTripletList(merged_local_data, tlist)
  CALL DestructMatrix(merged_local_data)

  !! Absolute Positions
  CALL ShiftTripletList(tlist, this%start_row - 1, this%start_column - 1)

  !! Figure out the length of the string for storing.
  triplet_list_string_length = 0
  DO II = 1, tlist%CurrentSize
#ifdef ISCOMPLEX
     CALL WriteMMLine(temp_string3, tlist%DATA(II)%index_row, &
          & tlist%DATA(II)%index_column, &
          & REAL(tlist%DATA(II)%point_value), &
          & AIMAG(tlist%DATA(II)%point_value), &
          & add_newline_in = .TRUE.)
#else
     CALL WriteMMLine(temp_string3, tlist%DATA(II)%index_row, &
          & tlist%DATA(II)%index_column, &
          & tlist%DATA(II)%point_value, add_newline_in = .TRUE.)
#endif
     WRITE(temp_string2, '(A)') ADJUSTL(temp_string3)
     triplet_list_string_length = triplet_list_string_length + &
          & LEN_TRIM(temp_string2)
     triplet_list_string_length = triplet_list_string_length + NEW_LINE_LENGTH
  END DO

  !! Write that string to the write buffer
  ALLOCATE(CHARACTER(LEN = triplet_list_string_length + 1) :: write_buffer)
  OFF_JJ = 1
  DO II = 1, tlist%CurrentSize
#ifdef ISCOMPLEX
     CALL WriteMMLine(temp_string3, tlist%DATA(II)%index_row, &
          & tlist%DATA(II)%index_column, REAL(tlist%DATA(II)%point_value), &
          & AIMAG(tlist%DATA(II)%point_value), add_newline_in = .TRUE.)
#else
     CALL WriteMMLine(temp_string3, tlist%DATA(II)%index_row, &
          & tlist%DATA(II)%index_column,tlist%DATA(II)%point_value, &
          & add_newline_in = .TRUE.)
#endif
     WRITE(temp_string2, '(A)') ADJUSTL(temp_string3)
     temp_length = LEN_TRIM(temp_string2) + NEW_LINE_LENGTH
     WRITE(write_buffer(OFF_JJ:OFF_JJ + temp_length), *) &
          & temp_string2(1:temp_length)
     OFF_JJ = OFF_JJ + temp_length
  END DO
  CALL DestructTripletList(tlist)

  !! Figure out the offset sizes
  ALLOCATE(local_values_buffer(this%process_grid%slice_size))
  CALL MPI_Allgather(triplet_list_string_length, 1, MPINTINTEGER,&
       & local_values_buffer, 1, MPINTINTEGER, &
       & this%process_grid%within_slice_comm, ierr)
  write_offset = 0
  write_offset = write_offset + header_size
  DO II = 1,this%process_grid%within_slice_rank
     write_offset = write_offset + local_values_buffer(II)
  END DO

  !! Global Write
  IF (this%process_grid%between_slice_rank .EQ. 0) THEN
     CALL MPI_File_open(this%process_grid%within_slice_comm, file_name, &
          & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY), MPI_INFO_NULL, &
          & mpi_file_handler,ierr)
     CALL MPI_File_set_size(mpi_file_handler, zero_size, ierr)

     !! Write Header
     header_offset = 0
     CALL MPI_File_set_view(mpi_file_handler, header_offset, MPI_CHARACTER, &
          & MPI_CHARACTER, "native", MPI_INFO_NULL, ierr)
     IF (this%process_grid%within_slice_rank .EQ. 0) THEN
        CALL MPI_File_write(mpi_file_handler, header_line1, &
             & LEN(header_line1), MPI_CHARACTER, message_status, ierr)
     END IF
     header_offset = header_offset + LEN(header_line1)
     CALL MPI_File_set_view(mpi_file_handler, header_offset, MPI_CHARACTER, &
          & MPI_CHARACTER, "native", MPI_INFO_NULL, ierr)
     IF (this%process_grid%within_slice_rank .EQ. 0) THEN
        CALL MPI_File_write(mpi_file_handler, header_line2, &
             & LEN(header_line2), MPI_CHARACTER, message_status, ierr)
     END IF

     !! Write Local Data
     CALL MPI_File_set_view(mpi_file_handler, write_offset, MPI_CHARACTER, &
          & MPI_CHARACTER, "native", MPI_INFO_NULL, ierr)
     CALL MPI_File_write(mpi_file_handler, write_buffer, &
          & triplet_list_string_length, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

     !! Cleanup
     CALL MPI_File_close(mpi_file_handler, ierr)
  END IF
  CALL MPI_Barrier(this%process_grid%global_comm, ierr)

  DEALLOCATE(header_line1)
  DEALLOCATE(header_line2)
  DEALLOCATE(write_buffer)
  DEALLOCATE(local_values_buffer)