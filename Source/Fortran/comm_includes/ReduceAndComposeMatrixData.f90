  !! Local Data
  INTEGER :: grid_error
  INTEGER :: II
  INTEGER :: total_values

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
  CALL MPI_IAllGatherv(matrix%inner_index,SIZE(matrix%values),MPI_INT, &
       & gathered_matrix%inner_index, helper%values_per_process, &
       & helper%displacement, MPI_INT, communicator, &
       & helper%inner_request, grid_error)
  CALL MPI_IAllGather(matrix%outer_index(2:), matrix%columns,&
       & MPI_INT, gathered_matrix%outer_index(2:), &
       & matrix%columns, MPI_INT, communicator, helper%outer_request, &
       & grid_error)
