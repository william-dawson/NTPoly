  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II

  CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

  !! Gather Information About Other Processes
  ALLOCATE(helper%values_per_process(helper%comm_size))
  ALLOCATE(helper%size_request_list(helper%comm_size*2))
  ALLOCATE(helper%size_status_list(helper%comm_size*2, MPI_STATUS_SIZE))
  DO II = 1, helper%comm_size
     CALL MPI_Isend(SIZE(matrix%values), 1, MPI_INT, II-1, 1, communicator, &
          & helper%size_request_list(II), grid_error)
     CALL MPI_Irecv(helper%values_per_process(II:II), 1, MPI_INT, II-1, 1, &
          & communicator, helper%size_request_list(helper%comm_size+II), &
          & grid_error)
  END DO
