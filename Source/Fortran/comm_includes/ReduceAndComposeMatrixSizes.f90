  !! Local Data
  INTEGER :: grid_error

  CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

  !! Build Storage
  CALL ConstructEmptyMatrix(gathered_matrix, &
       & matrix%rows,matrix%columns*helper%comm_size)
  gathered_matrix%outer_index(1) = 0

  !! Gather Information About Other Processes
  CALL MPI_IAllGather(matrix%outer_index(2:), matrix%columns,&
       & MPI_INT, gathered_matrix%outer_index(2:), &
       & matrix%columns, MPI_INT, communicator, helper%outer_request, &
       & grid_error)
