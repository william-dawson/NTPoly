  CALL ReduceAndSumMatrixSizes(matrix, comm, gathered_matrix, helper)
  DO WHILE(.NOT. TestReduceSizeRequest(helper))
  END DO

  CALL ReduceAndSumMatrixData(matrix, gathered_matrix, comm, helper)
  DO WHILE(.NOT. TestReduceInnerRequest(helper))
  END DO
  DO WHILE(.NOT. TestReduceDataRequest(helper))
  END DO
  
  CALL ReduceAndSumMatrixCleanup(matrix, gathered_matrix, threshold, helper)