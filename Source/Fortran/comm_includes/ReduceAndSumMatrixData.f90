  !! Local Data
  INTEGER :: grid_error
  INTEGER :: counter
  INTEGER :: sum_total_values, sum_outer_indices

  ALLOCATE(helper%displacement(helper%comm_size))

  !! Build Displacement List
  helper%displacement(1) = 0
  DO counter = 2, SIZE(helper%displacement)
     helper%displacement(counter) = helper%displacement(counter-1) + &
          & helper%values_per_process(counter-1)
  END DO

  !! Build Storage
  sum_total_values = SUM(helper%values_per_process)
  sum_outer_indices = (matrix%columns+1)*helper%comm_size
  CALL DestructMatrix(gathered_matrix)
  ALLOCATE(gathered_matrix%values(sum_total_values))
  ALLOCATE(gathered_matrix%inner_index(sum_total_values))
  ALLOCATE(gathered_matrix%outer_index(sum_outer_indices+1))

  !! MPI Calls
  CALL MPI_IAllGatherv(matrix%inner_index,SIZE(matrix%values),MPI_INT, &
       & gathered_matrix%inner_index, helper%values_per_process, &
       & helper%displacement, MPI_INT, communicator, helper%inner_request, &
       & grid_error)
  CALL MPI_IAllGather(matrix%outer_index, matrix%columns+1,&
       & MPI_INT, gathered_matrix%outer_index, matrix%columns+1, MPI_INT, &
       & communicator, helper%outer_request, grid_error)
