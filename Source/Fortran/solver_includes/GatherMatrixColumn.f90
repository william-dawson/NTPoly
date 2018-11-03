  !! Local Variables
  TYPE(ReduceHelper_t) :: gather_helper

  CALL TransposeMatrix(local_matrix, local_matrixT)
  CALL ReduceAndComposeMatrixSizes(local_matrixT, process_grid%column_comm, &
       & column_matrix, gather_helper)
  DO WHILE(.NOT. TestReduceSizeRequest(gather_helper))
  END DO
  CALL ReduceAndComposeMatrixData(local_matrixT, process_grid%column_comm, &
       & column_matrix, gather_helper)
  DO WHILE(.NOT. TestReduceInnerRequest(gather_helper))
  END DO
  DO WHILE(.NOT. TestReduceDataRequest(gather_helper))
  END DO
  CALL ReduceAndComposeMatrixCleanup(local_matrixT, column_matrix, &
       & gather_helper)

  CALL DestructMatrix(local_matrixT)
