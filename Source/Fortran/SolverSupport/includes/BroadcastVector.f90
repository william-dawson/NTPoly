  CALL MPI_Bcast(num_values, 1, MPI_INT, root, comm, err)
  CALL MPI_Bcast(indices(:num_values), num_values, MPI_INT, root, comm, err)
