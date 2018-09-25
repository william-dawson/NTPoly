  !! Local Data
  INTEGER :: grid_error
  INTEGER :: sum_outer_indices

  CALL MPI_Comm_size(communicator,helper%comm_size,grid_error)

  !! Build Storage
  CALL DestructMatrix(gathered_matrix)
  sum_outer_indices = (matrix%columns+1)*helper%comm_size
  ALLOCATE(gathered_matrix%outer_index(sum_outer_indices+1))

  !! Gather Outer Indices
  CALL MPI_IAllGather(matrix%outer_index, matrix%columns+1,&
       & MPI_INT, gathered_matrix%outer_index, matrix%columns+1, MPI_INT, &
       & communicator, helper%outer_request, grid_error)
