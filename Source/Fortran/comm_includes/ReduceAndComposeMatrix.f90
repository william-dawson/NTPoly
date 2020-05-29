  CALL ReduceAndComposeMatrixSizes(matrix, comm, gathered_matrix, helper)
  DO WHILE(.NOT. TestReduceSizeRequest(helper))
  END DO

  CALL ReduceAndComposeMatrixData(matrix, comm, gathered_matrix, helper)
  DO WHILE(.NOT. TestReduceInnerRequest(helper))
  END DO
  DO WHILE(.NOT. TestReduceDataRequest(helper))
  END DO

  CALL ReduceAndComposeMatrixCleanup(matrix, gathered_matrix, helper)