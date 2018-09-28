  CALL MPI_Bcast(num_values, 1, MPINTINTEGER, root, comm, err)
  CALL MPI_Bcast(indices(:num_values), num_values, MPINTINTEGER, root, &
       & comm, err)
