  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II
  INTEGER :: total_values
  INTEGER :: idx

  !! Compute values per process
  ALLOCATE(helper%values_per_process(helper%comm_size))
  DO II = 1, helper%comm_size
     idx = matrix%columns*II + 1
     helper%values_per_process(II) = gathered_matrix%outer_index(idx)
  END DO

  !! Build Displacement List
  ALLOCATE(helper%displacement(helper%comm_size))
  helper%displacement(1) = 0
  DO II = 2, SIZE(helper%displacement)
     helper%displacement(II) = helper%displacement(II-1) + &
          & helper%values_per_process(II-1)
  END DO

  !! Build Storage
  total_values = SUM(helper%values_per_process)
  ALLOCATE(gathered_matrix%values(total_values))
  ALLOCATE(gathered_matrix%inner_index(total_values))

  !! MPI Calls
  CALL MPI_IAllGatherv(matrix%inner_index,SIZE(matrix%values),MPINTINTEGER, &
       & gathered_matrix%inner_index, helper%values_per_process, &
       & helper%displacement, MPINTINTEGER, communicator, &
       & helper%inner_request, grid_error)
