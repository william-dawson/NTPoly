  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II
  INTEGER :: total_values
  INTEGER :: istart, iend, isize

  !! Send Receive Buffers
  ALLOCATE(helper%outer_request_list(helper%comm_size*2))
  ALLOCATE(helper%inner_request_list(helper%comm_size*2))
  ALLOCATE(helper%data_request_list(helper%comm_size*2))
  ALLOCATE(helper%outer_status_list(helper%comm_size*2, MPI_STATUS_SIZE))
  ALLOCATE(helper%inner_status_list(helper%comm_size*2, MPI_STATUS_SIZE))
  ALLOCATE(helper%data_status_list(helper%comm_size*2, MPI_STATUS_SIZE))

  !! Build Displacement List
  ALLOCATE(helper%displacement(helper%comm_size))
  helper%displacement(1) = 0
  DO II = 2, SIZE(helper%displacement)
     helper%displacement(II) = helper%displacement(II-1) + &
          & helper%values_per_process(II-1)
  END DO

  !! Build Storage
  CALL ConstructEmptyMatrix(gathered_matrix, &
       & matrix%rows,matrix%columns*helper%comm_size)
  total_values = SUM(helper%values_per_process)
  ALLOCATE(gathered_matrix%values(total_values))
  ALLOCATE(gathered_matrix%inner_index(total_values))
  gathered_matrix%outer_index(1) = 0

  !! MPI Calls
  DO II = 1, helper%comm_size
     !! Send/Recv inner index
     CALL MPI_ISend(matrix%inner_index, SIZE(matrix%inner_index), MPI_INT, &
          & II-1, 2, communicator, helper%inner_request_list(II), grid_error)
     istart = helper%displacement(II)+1
     isize = helper%values_per_process(II)
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%inner_index(istart:iend), isize, MPI_INT, &
          & II-1, 2, communicator, &
          & helper%inner_request_list(helper%comm_size+II), grid_error)
     !! Send/Recv Outer Index
     CALL MPI_ISend(matrix%outer_index(2:), matrix%columns, MPI_INT, &
          & II-1, 3, communicator, helper%outer_request_list(II), grid_error)
     istart = (matrix%columns)*(II-1)+2
     isize = matrix%columns
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%outer_index(istart:iend), isize, MPI_INT, &
          & II-1, 3, communicator, &
          & helper%outer_request_list(helper%comm_size+II), grid_error)
  END DO
