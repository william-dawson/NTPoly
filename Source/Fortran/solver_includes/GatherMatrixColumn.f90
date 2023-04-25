  CALL TransposeMatrix(local_matrix, local_matrixT)
  CALL ReduceAndComposeMatrix(local_matrixT, process_grid%column_comm, &
       & column_matrix)

  CALL DestructMatrix(local_matrixT)