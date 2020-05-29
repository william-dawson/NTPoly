  !! Local Data
  INTEGER, DIMENSION(:), ALLOCATABLE :: local_values_buffer
  INTEGER :: mpi_file_handler
  INTEGER(KIND=MPI_OFFSET_KIND) :: header_size
  INTEGER(KIND=MPI_OFFSET_KIND) :: write_offset
  !! Temporary Variables
  INTEGER :: bytes_per_int, bytes_per_entry
  INTEGER, DIMENSION(4) :: header_buffer
  INTEGER :: message_status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) :: zero_offset = 0
  INTEGER :: counter
  INTEGER :: ierr

  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data)

  !! Determine Write Location
  CALL MPI_Type_size(MPINTINTEGER, bytes_per_int, ierr)
  CALL MPI_Type_extent(triplet_mpi_type, bytes_per_entry, ierr)
  header_size = bytes_per_int*4
  ALLOCATE(local_values_buffer(this%process_grid%slice_size))
  CALL MPI_Allgather(SIZE(merged_local_data%values), 1, MPINTINTEGER,&
       & local_values_buffer, 1, MPINTINTEGER,&
       & this%process_grid%within_slice_comm,ierr)
  write_offset = 0
  write_offset = write_offset + header_size
  DO counter = 1,this%process_grid%within_slice_rank
     write_offset = write_offset + &
          & local_values_buffer(counter)*(bytes_per_entry)
  END DO

  !! Write The File
  IF (this%process_grid%between_slice_rank .EQ. 0) THEN
     !! Create Special MPI Type
     CALL MatrixToTripletList(merged_local_data, triplet_list)
     !! Absolute Positions
     CALL ShiftTripletList(triplet_list, this%start_row - 1, &
          & this%start_column - 1)
     CALL MPI_File_open(this%process_grid%within_slice_comm, file_name,&
          & IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY), MPI_INFO_NULL, &
          & mpi_file_handler, ierr)
     !! Write Header
     IF (this%process_grid%within_slice_rank .EQ. 0) THEN
        header_buffer(1) = this%actual_matrix_dimension
        header_buffer(2) = this%actual_matrix_dimension
        header_buffer(3) = SUM(local_values_buffer)
        IF (this%is_complex) THEN
           header_buffer(4) = 1
        ELSE
           header_buffer(4) = 0
        END IF
        CALL MPI_File_write_at(mpi_file_handler, zero_offset, header_buffer, &
             & 4, MPINTINTEGER, message_status, ierr)
     END IF
     !! Write The Rest
     CALL MPI_File_set_view(mpi_file_handler,write_offset,triplet_mpi_type,&
          & triplet_mpi_type,"native",MPI_INFO_NULL,ierr)
     CALL MPI_File_write(mpi_file_handler, triplet_list%DATA, &
          & triplet_list%CurrentSize, triplet_mpi_type, MPI_STATUS_IGNORE, &
          & ierr)

     !! Cleanup
     CALL MPI_File_close(mpi_file_handler,ierr)
     CALL DestructTripletList(triplet_list)
  END IF
  DEALLOCATE(local_values_buffer)
  CALL MPI_Barrier(this%process_grid%global_comm,ierr)
  CALL DestructMatrix(merged_local_data)
