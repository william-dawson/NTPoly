  !! Local Variables
  TYPE(ReduceHelper_t) :: gather_helper
  INTEGER :: mpi_status(MPI_STATUS_SIZE)
  INTEGER :: mpi_err

  CALL TransposeMatrix(local_matrix, local_matrixT)
  CALL ReduceAndComposeMatrixSizes(local_matrixT, process_grid%column_comm, &
       & column_matrix, gather_helper)
  CALL MPI_Wait(gather_helper%outer_request, mpi_status, mpi_err)
  CALL ReduceAndComposeMatrixData(local_matrixT, process_grid%column_comm, &
       & column_matrix, gather_helper)
  CALL MPI_Wait(gather_helper%inner_request,mpi_status,mpi_err)
  CALL MPI_Wait(gather_helper%data_request,mpi_status,mpi_err)
  CALL ReduceAndComposeMatrixCleanup(local_matrixT, column_matrix, &
       & gather_helper)

  CALL DestructMatrix(local_matrixT)
