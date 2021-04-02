  CALL TransposeMatrix(local_matrix, local_matrixT)
  CALL ReduceAndComposeMatrix(local_matrixT, column_matrix, &
       & process_grid%column_comm)

  CALL DestructMatrix(local_matrixT)