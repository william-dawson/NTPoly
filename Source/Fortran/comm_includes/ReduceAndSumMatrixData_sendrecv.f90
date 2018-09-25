  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II
  INTEGER :: sum_total_values, sum_outer_indices
  INTEGER :: istart, isize, iend

  !! Send Receive Buffers
  ALLOCATE(helper%outer_send_request_list(helper%comm_size))
  ALLOCATE(helper%outer_recv_request_list(helper%comm_size))
  ALLOCATE(helper%inner_send_request_list(helper%comm_size))
  ALLOCATE(helper%inner_recv_request_list(helper%comm_size))
  ALLOCATE(helper%data_send_request_list(helper%comm_size))
  ALLOCATE(helper%data_recv_request_list(helper%comm_size))

  ALLOCATE(helper%displacement(helper%comm_size))
  !! Build Displacement List
  helper%displacement(1) = 0
  DO II = 2, SIZE(helper%displacement)
     helper%displacement(II) = helper%displacement(II-1) + &
          & helper%values_per_process(II-1)
  END DO

  !! Build Storage
  sum_total_values = SUM(helper%values_per_process)
  sum_outer_indices = (matrix%columns+1)*helper%comm_size
  CALL DestructMatrix(gathered_matrix)
  ALLOCATE(gathered_matrix%values(sum_total_values))
  ALLOCATE(gathered_matrix%inner_index(sum_total_values))
  ALLOCATE(gathered_matrix%outer_index(sum_outer_indices+1))

  !! MPI Calls
  DO II = 1, helper%comm_size
     !! Send/Recv inner index
     CALL MPI_ISend(matrix%inner_index, SIZE(matrix%inner_index), MPI_INT, &
          & II-1, 2, communicator, helper%inner_send_request_list(II), &
          & grid_error)
     istart = helper%displacement(II)+1
     isize = helper%values_per_process(II)
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%inner_index(istart:iend), isize, MPI_INT, &
          & II-1, 2, communicator, &
          & helper%inner_recv_request_list(II), grid_error)
  END DO
