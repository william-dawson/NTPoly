  !! Local Data
  INTEGER :: sum_outer_indices
  INTEGER :: ierr

  CALL MPI_Comm_size(comm, helper%comm_size, ierr)

  !! Build Storage
  CALL DestructMatrix(gathered_matrix)
  sum_outer_indices = (matrix%columns + 1) * helper%comm_size
  ALLOCATE(gathered_matrix%outer_index(sum_outer_indices + 1))

  !! Gather Outer Indices
  CALL MPI_IAllGather(matrix%outer_index, matrix%columns+1, &
       & MPINTINTEGER, gathered_matrix%outer_index, matrix%columns + 1, &
       & MPINTINTEGER, comm, helper%outer_request, ierr)
