  !! Local Data
  INTEGER :: ierr

  CALL MPI_Comm_size(comm, helper%comm_size, ierr)

  !! Build Storage
  CALL ConstructEmptyMatrix(gathered_matrix, &
       & matrix%rows, matrix%columns * helper%comm_size)
  gathered_matrix%outer_index(1) = 0

  !! Gather Information About Other Processes
  CALL MPI_IAllGather(matrix%outer_index(2:), matrix%columns,&
       & MPINTINTEGER, gathered_matrix%outer_index(2:), &
       & matrix%columns, MPINTINTEGER, comm, helper%outer_request, ierr)
