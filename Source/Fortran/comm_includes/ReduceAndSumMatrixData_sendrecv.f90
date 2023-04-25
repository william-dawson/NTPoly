  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II
  INTEGER :: sum_total_values
  INTEGER :: istart, isize, iend
  INTEGER :: idx

  !! Send Receive Buffers
  ALLOCATE(helper%inner_send_request_list(helper%comm_size))
  ALLOCATE(helper%inner_recv_request_list(helper%comm_size))
  ALLOCATE(helper%data_send_request_list(helper%comm_size))
  ALLOCATE(helper%data_recv_request_list(helper%comm_size))

  !! Compute values per process
  ALLOCATE(helper%values_per_process(helper%comm_size))
  DO II = 1, helper%comm_size
     idx = (matrix%columns + 1) * II
     helper%values_per_process(II) = gathered_matrix%outer_index(idx)
  END DO

  !! Build Displacement List
  ALLOCATE(helper%displacement(helper%comm_size))
  helper%displacement(1) = 0
  DO II = 2, SIZE(helper%displacement)
     helper%displacement(II) = helper%displacement(II - 1) + &
          & helper%values_per_process(II - 1)
  END DO

  !! Build Storage
  sum_total_values = SUM(helper%values_per_process)
  ALLOCATE(gathered_matrix%values(sum_total_values))
  ALLOCATE(gathered_matrix%inner_index(sum_total_values))

  !! MPI Calls
  DO II = 1, helper%comm_size
     !! Send/Recv inner index
     CALL MPI_ISend(matrix%inner_index, SIZE(matrix%inner_index), &
          & MPINTINTEGER, II - 1, 2, communicator, &
          & helper%inner_send_request_list(II), grid_error)
     istart = helper%displacement(II) + 1
     isize = helper%values_per_process(II)
     iend = istart + isize - 1
     CALL MPI_Irecv(gathered_matrix%inner_index(istart:iend), isize, &
          & MPINTINTEGER, II - 1, 2, communicator, &
          & helper%inner_recv_request_list(II), grid_error)
  END DO
