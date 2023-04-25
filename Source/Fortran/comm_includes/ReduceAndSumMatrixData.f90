  !! Local Data
  INTEGER :: II, idx
  INTEGER :: sum_total_values
  INTEGER :: ierr

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
  CALL MPI_IAllGatherv(matrix%inner_index, SIZE(matrix%values), MPINTINTEGER, &
       & gathered_matrix%inner_index, helper%values_per_process, &
       & helper%displacement, MPINTINTEGER, comm, helper%inner_request, ierr)
